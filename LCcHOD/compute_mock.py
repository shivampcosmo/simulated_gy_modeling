#!/usr/bin/env python
from os import path
import argparse
import subprocess
import numpy as np
import scipy.optimize
import h5py as h5

import config
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
    
    #print(logMmin)
    #print(np.sum(ngal*dM))
    return np.sum(ngal*dM)

def compute_HOD_parameters(ndens=None, nfid=None, M1_over_Mmin=None, M0_over_M1=None, alpha=None, siglogM=None, f_cen=None, halos=None, header=None, logMmin_guess=12.5):
    """
    Compute the physical (er, Msun h^-1) HOD parameters from ratios of masses
    and the desired number density of galaxies.
    We convert the number density of galaxies into Mmin using the halo mass function.
    """    
    #logM = np.linspace(10.0, 16.0, 512)
    logM = np.linspace(12.0, 16.0, 512) #this is temporary to investigate halo model calculations
    
    with h5.File(halos, mode='r') as catalog:
        bin_counts, bin_edges = np.histogram(np.log10(catalog['halos']['mass']),bins=logM)

    cf = config.AbacusConfigFile(header)
    vol = cf.boxSize**3
    dM = np.diff(10**(logM))
    mass_function = bin_counts / vol / dM
    
    logM = logM[:-1]
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
                     omega_m,redshift,boxsize,siglogM,logMmin,logM0,logM1,alpha,q_cen,q_sat,A_con,f_cen,sigphot,seed=None):
    script_path = path.dirname(path.abspath(__file__))+"/../cHOD/compute_mocks"
    cmd_line = [script_path,str(omega_m),str(redshift),\
                    str(siglogM),str(logMmin),str(logM0),str(logM1),\
                    str(alpha),str(q_cen),str(q_sat),str(A_con),str(f_cen),str(sigphot),str(boxsize),str(seed),\
                    halo_file,galaxy_file,env_file]
    subprocess.call(cmd_line)

def compute(halo_file,env_file,header_file,output_file,siglogM=None,logMmin=None,logM0=None,logM1=None,alpha=None,q_cen=None,q_sat=None,A_con=None,f_cen=None,sigphot=None,seed=42,njackknife=8):
    cf = config.AbacusConfigFile(header_file)
    boxsize = cf.boxSize
    omeganow_m = cf.OmegaNow_m
    omega_m = cf.Omega_M
    redshift = cf.redshift
    
    ## now compute HOD
    populate_hod(halo_file, output_file, env_file, omega_m,redshift,boxsize,\
                         siglogM,logMmin,logM0,logM1,alpha,q_cen,q_sat,A_con,f_cen,sigphot,seed=seed)

parser = argparse.ArgumentParser()
parser.add_argument('hod_params_path')
parser.add_argument('header_path')
parser.add_argument('halo_path')
parser.add_argument('output_path')
parser.add_argument('env_path')

args = parser.parse_args()

# read meta-HOD parameters
myconfigparser = configparser.ConfigParser()
myconfigparser.read(args.hod_params_path)
params = myconfigparser['params']

if 'fracngal' in params.keys():
    fracngal = float(params['fracngal'])
else:
    fracngal = 1.000
    
# find HOD parameters
logMmin, logM0, logM1 = compute_HOD_parameters(ndens=fracngal,nfid=float(params['ngal']),siglogM=float(params['siglogm']),M0_over_M1=float(params['m0_over_m1']),M1_over_Mmin=float(params['m1_over_mmin']),alpha=float(params['alpha']),f_cen=float(params['f_cen']),halos=args.halo_path,header=args.header_path)

print(fracngal * float(params['ngal']))
print(logMmin)
print(logM0)
print(logM1)

#output_prefix = "HOD_%.2f_%.2f_%.2f_%.2f_%.2f_seed_%s" % (siglogM,logMmin,logM0,logM1,alpha,seed)
#output_prefix = "NHOD_%f_%f_%f_%f_%f_seed_%s" % (args.siglogM,args.ngal,args.M0_over_M1,args.M1_over_Mmin,args.alpha,seed)
seed = int(params['seed'])

import sys
print("\tq_cen: ",params['q_cen'], file=sys.stderr)

compute(args.halo_path,args.env_path,args.header_path,args.output_path,
        params['siglogm'],logMmin,logM0,logM1,params['alpha'],params['q_cen'],params['q_sat'],params['A_con'],params['f_cen'],params['sigphot'],
        seed=seed)

