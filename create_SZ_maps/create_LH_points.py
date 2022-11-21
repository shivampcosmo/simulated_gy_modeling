import numpy as np
import matplotlib.pyplot as plt

from smt.sampling_methods import LHS

import copy
import numpy as np
import configparser
# from pyDOE import lhs
config = configparser.ConfigParser()

fname = 'params_all.ini'
config.read(fname)


sec_names = list(config.keys())[1:]

all_vary_names = []
all_vary_minv = []
all_vary_maxv = []
for sec in sec_names:
    var_names = list(config[sec].keys())
    for var in var_names:
        rangev = (config[sec][var])
        rangev_arr = list(map(float,rangev.split())) 
        if len(rangev_arr) > 1:
            all_vary_names.append(sec + '--' + var)
            all_vary_minv.append(rangev_arr[0])
            all_vary_maxv.append(rangev_arr[2])        

    
nvar_all = len(all_vary_names)
xlimits = np.zeros((nvar_all,2))
for jv in range((nvar_all)):
    xlimits[jv,0] = 0.0
    xlimits[jv,1] = 1.0

ntot = 2
rand_states = np.linspace(104,1004,ntot).astype(int)
rsv = np.arange(0,ntot).astype(int)

for jsv in range(len(rsv)):

    # sampling = LHS(xlimits=xlimits,criterion='ese',random_state=int(rand_states[jsv]))
    sampling = LHS(xlimits=xlimits,criterion='cm',random_state=int(rand_states[jsv]))

    num = 100
    x = sampling(num)

    nLHS = num
    LHS_points = x

    LHS_ids = np.arange(nLHS).astype(int)

    LHS_points_final = np.zeros_like(LHS_points)
    for jv in range(nvar_all):
        LHS_points_final[:,jv] = all_vary_minv[jv] + (all_vary_maxv[jv]-all_vary_minv[jv]) * LHS_points[:,jv]

    LHS_points_final = np.hstack((np.array([LHS_ids]).T,LHS_points_final))

    first_line = 'emusave--lhsid    '
    for var in all_vary_names:
        first_line += str(var) + '    '

    np.savetxt('sample_chain_vary_params_np_' + str(nvar_all) + '_rs' + str(rsv[jsv]) + '.txt',LHS_points_final,header=first_line)

