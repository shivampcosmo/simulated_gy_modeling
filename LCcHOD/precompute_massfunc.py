import argparse
import numpy as np
import h5py as h5
import astropy.cosmology

def compute_mass_function(zmin, zmax, infile, V):
    infile = h5.File(infile, 'r')

    halos = infile['halos']

    print(len(halos))

    all_masses = halos['M_interp']
    redshift_interp = halos['redshift_interp']

    infile.close()

    index = np.where( (redshift_interp > zmin) & (redshift_interp < zmax) )

    print(len(index))

    binned_masses = all_masses[index]

    mass_bins = np.logspace(11, 15, 101)
    mass_avg = (mass_bins[1:] + mass_bins[:-1])/2.0

    hist, bin_edges = np.histogram(binned_masses, mass_bins)

    dM = mass_bins[1:] - mass_bins[:-1]

    dndm = hist / (dM * V)
    
    return mass_avg, dM, dndm

parser = argparse.ArgumentParser()
parser.add_argument('halo_file')
parser.add_argument('zmin')
parser.add_argument('zmax')
parser.add_argument('outfile')
args = parser.parse_args()

cosmo = astropy.cosmology.FlatLambdaCDM(H0 = 67.36, Om0 = 0.314, Neff = 3.04) #currently hardcoded to Summit cosmology

zmin = float(args.zmin)
zmax = float(args.zmax)

Rmin = cosmo.comoving_distance(zmin).value * cosmo.h #comoving distance in Mpc/h
Rmax = cosmo.comoving_distance(zmax).value * cosmo.h 

V = (4. / 3.) * np.pi * (Rmax**3.0 - Rmin**3.0) * (1./8.) #hardcoded for Summit lightcone geometry, full octant to z ~ 0.8

M, dM, dndM = compute_mass_function(zmin, zmax, args.halo_file, V)

np.savetxt(args.outfile, np.transpose(np.array([M, dM, dndM])))
