import os, sys
import numpy as np
import matplotlib.pyplot as plt

def plot_all_clones_combined(proj_dir, sys, Clones_list, swaps_per_ns, note, verbose=True):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4), sharex=True, constrained_layout=True)
    cmap = plt.get_cmap('tab20')

    # ========== Subplot 1: Free Energy Trajectories ==========
    for c in Clones_list:
        color = cmap(c % 20)
        energy_path = f'{proj_dir}/EE_analysis/scaped_data/{sys}_CLONE{c}_feb.npy'
        if not os.path.exists(energy_path):
            if verbose:
                print(f"File not found: {energy_path}")
            continue
        clone_frame_energies = np.load(energy_path)
        ax1.scatter(np.arange(len(clone_frame_energies)) / swaps_per_ns,
                    clone_frame_energies[:, -1], s=1, color=color, label=f"CLONE {c}")
    ax1.set_ylabel(r'$\Delta G_{\mathrm{L}}$ (kT)')
    ax1.set_ylim(30, 50) #for water
    ax1.legend(loc='lower right', ncol=2, fontsize='5')
    #ax1.set_ylim(40, 50) #for complex
    #ax1.set_ylim(-5, 7) #for A16F
    ax1.set_xlim(-1, 30)
    #ax1.set_title('Free Energy Trajectories (All Clones)')

    # ========== Subplot 2: WL-Increment ==========
    for c in Clones_list:
        increment_path = f'{proj_dir}/EE_analysis/scaped_data/{sys}_CLONE{c}_all_increments.npy'
        if not os.path.exists(increment_path):
            if verbose:
                print(f"File not found: {increment_path}")
            continue
        all_increments = np.load(increment_path)
        log_increments = np.log10(all_increments.astype(np.float64))
        color = cmap(c % 20)
        ax2.scatter(np.arange(len(log_increments)) / swaps_per_ns,
                    log_increments, s=1, color=color, label=f"CLONE {c}")
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel(r'log$_{10}$ WL Increment')
    ax2.set_ylim(-4, 1)
    ax2.set_xlim(-1, 200)
    #ax2.set_xlim(-5, 420)
    #ax2.set_title('WL Increment Trajectories (All Clones)')

    # ========== Save the Figure ==========
    output_path = f'{proj_dir}/EE_analysis/plots/{sys}_{note}_combined_EE_WL_plot.png'
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    if verbose:
        print(f"✅ Plot saved to {output_path}")

sys.path.append('/home/tun50867/work/git/codes/system_preparation')
from scape_log_plot_dG_calEE import *

proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45'
#proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/A16F_relative'

#sys = 'Simulation_Files'
#Clones_list = [1,2,3,4,5,6] #for complex 1/t
#Clones_list = [7,8,9,10,11] #for complex no 1/t

sys = 'water'
#Clones_list = [0,1,3,4] #for water 1/t
#Clones_list = [2,5,7,8,9] #for water no 1/t
#Clones_list = [10,11,12,13,14,15,16,17,18,19]
Clones_list = [20,21,22,23,24,25,26,27,28,29]
#Clones_list = [40,41,42,43,44,45,46,47,48,49]
#sys = 'complex'
#Clones_list = [5,6,7,8,9] #for complex no 1/t
#Clones_list = [0,1,2,3,4] #for complex 1/t
mdp_file = f'{proj_dir}/{sys}/CLONE{Clones_list[0]}/frame1.mdp'
swaps_per_ns, lambdas, eq_dist, k = scrape_mdp(mdp_file)

plot_all_clones_combined(proj_dir, sys, Clones_list, swaps_per_ns,note = 'C20_29_1overt', verbose=True)

