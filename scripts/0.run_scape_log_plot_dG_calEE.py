
import glob, os, collections, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec, rc

sys.path.append('/home/tun50867/work/git/codes/system_preparation')
from scape_log_plot_dG_calEE import *

#proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45'
#proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/A16F_relative'
proj_dir = '/home/tun50867/work/git/one_over_t_project/systems'

#RUNs = ['Simulation_Files']
#RUNs = ['water']
#RUNs = ['complex']

RUNs = ['Paracetamol']

all_RUNs_dG =[]
for sys in RUNs : # for nRUNs
    print(f'for sys_name: {sys}')
    Clones_list = [0,1,2,3,4,5,6,7,8]
    #Clones_list = [20,21,22,23,24,25,26,27,28,29]
    #Clones_list = [40,41,42,43,44,45,46,47,48,49]
    for c in Clones_list: #for nCLONEs
        mdp_file = f'{proj_dir}/{sys}/CLONE{c}/frame1.mdp'
        swaps_per_ns, lambdas, eq_dist, k = scrape_mdp(mdp_file)

        nlog_files = 15
        if nlog_files != 0: 
            scrape_log(proj_dir, sys, c, log_list = range(0,nlog_files),stride=10, verbose=True)
        plot(proj_dir, sys, c, swaps_per_ns, verbose=True)
    
    #plot_all_clones(proj_dir, sys, Clones_list, swaps_per_ns,note = '1overt', verbose=True)

    #weighted_mean_avg_G = cal_weighted_mean_avg_dG(proj_dir, sys, Clones_list)
    #all_RUNs_dG.append(weighted_mean_avg_G)
#np.save(f'{proj_dir}/EE_analysis/plots/all_syss_dG.npy', all_RUNs_dG)
