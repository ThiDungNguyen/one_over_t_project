import sys
sys.path.append('/home/tun50867/work/git/codes/system_preparation')
from scape_log_plot_dG_calEE import *
import matplotlib.pyplot as plt
# === Apply scientific style settings ===
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.linewidth": 1.2,
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "legend.fontsize": 10,
    "lines.linewidth": 2,
})

# === Load Data ===
def load_clones(clone_list):
    free_energies = []
    swaps_per_ns_list = []
    for c in clone_list:
        swaps_per_ns, lambdas, eq_dist, k = scrape_mdp(mdp_file = f'{proj_dir}/{system_name}/CLONE{c}/frame1.mdp')
        swaps_per_ns_list.append(swaps_per_ns)
        path = f'{proj_dir}/EE_analysis/scaped_data/{system_name}_CLONE{c}_feb.npy'
        if os.path.exists(path):
            data = np.load(path)
            free_energies.append(data[:, -1])  # last column: free energy at each time
        else:
            print(f"File not found: {path}")
        #print('swaps_per_ns_list',swaps_per_ns_list)
    return free_energies, swaps_per_ns_list[0]


def plot_free_energy_comparison(free_energies_1t, free_energies_no1t, time_ns, label_1t="1/t tempering", label_no1t="No 1/t", color_1g='red', color_no1t='black'):
    # Convert to arrays
    fe_1t = np.array(free_energies_1t)     # Shape: (n_clones, n_timepoints)
    fe_no1t = np.array(free_energies_no1t)

    # Mean and std across clones
    mean_1t = np.mean(fe_1t, axis=0)
    std_1t = np.std(fe_1t, axis=0)

    mean_no1t = np.mean(fe_no1t, axis=0)
    std_no1t = np.std(fe_no1t, axis=0)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, mean_1t, color='C0', label="1/t tempering")
    plt.fill_between(time_ns, mean_1t - std_1t, mean_1t + std_1t, color="C0", alpha=0.3)

    plt.plot(time_ns, mean_no1t, color='C1', label="No 1/t")
    plt.fill_between(time_ns, mean_no1t - std_no1t, mean_no1t + std_no1t, color="C1", alpha=0.3)

    # Final point with error bar
    plt.errorbar(time_ns[-1], mean_1t[-1], yerr=std_1t[-1], fmt='o', color='C0', capsize=5)
    plt.errorbar(time_ns[-1], mean_no1t[-1], yerr=std_no1t[-1], fmt='o', color='C1', capsize=5)

    # Labels
    plt.xlabel("Simulation Time (ns)", fontsize=14)
    plt.ylabel("Free Energy (ΔG)", fontsize=14)
    #plt.title("Comparison of ΔG Over Time: 1/t vs No 1/t", fontsize=14)
    plt.legend(loc='upper right')
    plt.ylim(38, 52)
    plt.tight_layout()
    plt.savefig(f'{proj_dir}/EE_analysis/plots/{system_name}_average_EE_result.pdf', format='pdf', dpi=600, bbox_inches='tight')
    plt.close()



#proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45'
#system_name = 'complex'
#list_clones_1t = [0,1,2,3,4]
#list_clones_no1t = [5,6,7,8,9]
#system_name = 'water'
#list_clones_1t = [0,1,3,4]
#list_clones_no1t = [2,5]

proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/Paracetamol'
system_name = 'Paracetamol'
list_clones_no1t = [0,1,2,3,4,5,6,7,8,9]

#swaps_per_ns, lambdas, eq_dist, k = scrape_mdp(mdp_file = f'{proj_dir}/{sys}/CLONE{c}/frame1.mdp')

#free_energies_1t, swaps_per_ns_1t = load_clones(list_clones_1t)
free_energies_no1t, swaps_per_ns_no1t = load_clones(list_clones_no1t)

swaps_per_ns = swaps_per_ns_no1t  # safe to use now
n_steps = len(free_energies_1t[0])
time_ns = np.arange(n_steps) / swaps_per_ns

# Align all trajectories to the shortest length
min_length = min([len(f) for f in free_energies_1t + free_energies_no1t])

#free_energies_1t = [f[:min_length] for f in free_energies_1t]
free_energies_no1t = [f[:min_length] for f in free_energies_no1t]
time_ns = time_ns[:min_length]

## Step 1: Find index for 20 ns
#index_5ns = int(5 * swaps_per_ns)
#
## Step 2: Truncate all data to after 20 ns
#fe_1t = [f[index_5ns:] for f in free_energies_1t]
#fe_no1t = [f[index_5ns:] for f in free_energies_no1t]
#time_ns_truncated = time_ns[index_5ns:]
#
## Step 3: Align lengths
#min_length = min([len(f) for f in fe_1t + fe_no1t])
#fe_1t = [f[:min_length] for f in fe_1t]
#fe_no1t = [f[:min_length] for f in fe_no1t]
#time_ns_truncated = time_ns_truncated[:min_length]
#
## Step 4: Plot
#plot_free_energy_comparison(fe_1t, fe_no1t, time_ns_truncated)

plot_free_energy_comparison(free_energies_1t, free_energies_no1t, time_ns)

