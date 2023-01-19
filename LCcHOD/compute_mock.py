#!/usr/bin/env python
from os import path
import argparse
import subprocess
import numpy as np
import scipy.optimize
import h5py as h5

import configparser

def compute_ngal(logM, dM, mass_function, logMmin, M1_over_Mmin, M0_over_M1, alpha, siglogM, f_cen):
    M = 10.**logM
    Mmin = 10.**logMmin
    M1 = M1_over_Mmin*Mmin
    M0 = M0_over_M1*M1
    mean_cen = f_cen * 0.5 * (1. + scipy.special.erf((logM - logMmin)/siglogM))
    mean_sat = 0.5 * (1. + scipy.special.erf((logM - logMmin)/siglogM)) * (((M - M0) / M1)**alpha)

    ngal = np.zeros(logM.shape)
    ngal += mean_cen
    ngal[logM > np.log10(M0)] += mean_sat[logM > np.log10(M0)]
    ngal *= mass_function

    return np.sum(ngal*dM)

def compute_HOD_parameters(ndens=None, nfid=None, M1_over_Mmin=None, M0_over_M1=None, alpha=None, siglogM=None, f_cen=None, dn_dm_file=None, logMmin_guess=12.5):
    """
    Compute the physical (er, Msun h^-1) HOD parameters from ratios of masses
    and the desired number density of galaxies.
    We convert the number density of galaxies into Mmin using the halo mass function.
    """    
    
    logM, dM, mass_function = np.genfromtxt(dn_dm_file)

    """
    now solve for logMmin:
    compute_ngal(...) - ndens*nfid = 0
    """
    objective = lambda logMmin: compute_ngal(logM, dM, mass_function, logMmin, M1_over_Mmin, M0_over_M1, alpha, siglogM, f_cen) - ndens*nfid
    
    logMmin = scipy.optimize.newton(objective, logMmin_guess, maxiter = 1000)
    Mmin = 10**(logMmin)
    M1 = Mmin*M1_over_Mmin
    M0 = M1*M0_over_M1
    
    return logMmin, np.log10(M0), np.log10(M1)

def populate_hod(halo_file, galaxy_file, env_file, \
                 siglogM,logMmin,logM0,logM1,alpha,q_cen,q_sat,A_con,f_cen,z_pivot,slope,seed=None):
    script_path = path.dirname(path.abspath(__file__))+"/compute_mocks"
    cmd_line = [script_path,\
                    str(siglogM),str(logMmin),str(logM0),str(logM1),\
                    str(alpha),str(f_cen),str(q_cen),str(q_sat),str(A_con),str(z_pivot),str(slope),str(seed),\
                    halo_file,galaxy_file,env_file]
    subprocess.call(cmd_line)

parser = argparse.ArgumentParser()
parser.add_argument('hod_param_LHS_table')
parser.add_argument('param_index')
parser.add_argument('precomp_massfunc')
parser.add_argument('halo_path')
parser.add_argument('output_path')
parser.add_argument('env_path')
args = parser.parse_args()

#Need to fill in: open parameter table for Latin-Hypercube Sample, read in the requisite row and read in the hod parameters.
#fracngal, nfid, siglogM, M1oMmin, M0oM1, alpha, fcen, Qcen, Qsat, Acon, z_pivot, slope ... etc.
#We may fix some of these parameters.

# find HOD parameters
logMmin, logM0, logM1 = compute_HOD_parameters(ndens=fracngal,nfid=nfid,siglogM=siglogM,M0_over_M1=M0oM1,M1_over_Mmin=M1oMmin,alpha=alpha,f_cen=fcen,dn_dm_file=args.precomp_massfunc)

print(fracngal * nfid)
print(logMmin)
print(logM0)
print(logM1)

seed = 42 #use non-zero seed to be safe

populate_hod(args.halo_path,args.output_path,args.env_path,
        siglogM,logMmin,logM0,logM1,alpha,Qcen,Qsat,Acon,fcen,z_pivot,slope,
        seed=seed)

