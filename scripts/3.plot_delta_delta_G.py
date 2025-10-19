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

def cal_weighted_mean_avg_dG(proj_dir, sys, Clones_list):
    all_clones_dG_var, all_clones_dG = [] , []
    for c in Clones_list:
        if not os.path.exists(f'{proj_dir}/EE_analysis/scaped_data/{sys}_CLONE{c}_feb.npy'):
            dG, dG_sigma, dG_var = np.nan, np.nan, np.nan
        else:
            clone_frame_energies = np.load(f'{proj_dir}/EE_analysis/scaped_data/{sys}_CLONE{c}_feb.npy')[:,-1]
            ihalf = int(len(clone_frame_energies)/2)
            dG, dG_sigma, dG_var = np.average(clone_frame_energies[ihalf:]), np.std(clone_frame_energies[ihalf:])+1e-5, np.var(clone_frame_energies[ihalf:])+1e-5 # to avoid zeros
        all_clones_dG_var.append(dG_var)
        all_clones_dG.append(dG)
        print('dG, dG_sigma, dG_var',dG, dG_sigma, dG_var)

    # compute 1/var^2-weighted mean G
    np.save(f'{proj_dir}/EE_analysis/plots/{sys}_all_clones_dG_var.npy', all_clones_dG_var)
    np.save(f'{proj_dir}/EE_analysis/plots/{sys}_all_clones_dG.npy', all_clones_dG)

    # Filter out NaN values
    valid_dG_var = np.array([var for var in all_clones_dG_var if not np.isnan(var)])
    valid_dG = np.array([dG for dG in all_clones_dG if not np.isnan(dG)])

    # Compute the weighted mean average of dG
    if len(valid_dG_var) == 0 or len(valid_dG) == 0:
        weighted_mean_avg_G = np.nan
    else:
        inv_var_G = 1. / valid_dG_var
        weighted_mean_avg_G = np.sum(inv_var_G * valid_dG) / np.sum(inv_var_G)
    #print('inv_var_G',inv_var_G)
    #print('weighted_mean_avg_G',weighted_mean_avg_G)
    return weighted_mean_avg_G


def cal_per_clone_avg_dG(proj_dir, sys, Clones_list, half_data=True, save=True, verbose=True):
    all_clones_dG = []
    all_clones_dG_var = []

    for c in Clones_list:
        npy_path = f'{proj_dir}/EE_analysis/scaped_data/{sys}_CLONE{c}_feb.npy'
        if not os.path.exists(npy_path):
            dG, dG_sigma, dG_var = np.nan, np.nan, np.nan
        else:
            clone_frame_energies = np.load(npy_path)[:, -1]
            data_to_use = clone_frame_energies[int(len(clone_frame_energies)/2):] if half_data else clone_frame_energies

            dG = np.mean(data_to_use)
            dG_sigma = np.std(data_to_use) + 1e-5
            dG_var = np.var(data_to_use) + 1e-5

        all_clones_dG.append(dG)
        all_clones_dG_var.append(dG_var)

        if verbose:
            print(f"CLONE {c}: dG = {dG:.4f}, σ = {dG_sigma:.4f}, var = {dG_var:.4f}")

    if save:
        plot_dir = f'{proj_dir}/EE_analysis/plots'
        os.makedirs(plot_dir, exist_ok=True)
        np.save(f'{plot_dir}/{sys}_all_clones_dG.npy', all_clones_dG)
        np.save(f'{plot_dir}/{sys}_all_clones_dG_var.npy', all_clones_dG_var)

    return all_clones_dG, all_clones_dG_var

def compute_weighted_mean_dG(dG_list, var_list):
    dG = np.array(dG_list)
    var = np.array(var_list)

    valid = ~np.isnan(dG) & ~np.isnan(var)
    if not np.any(valid):
        return np.nan

    inv_var = 1.0 / var[valid]
    weighted_mean = np.sum(inv_var * dG[valid]) / np.sum(inv_var)
    return weighted_mean

def compute_weighted_average(dG_list, var_list):
    dG = np.array(dG_list)
    var = np.array(var_list)

    # Filter out NaNs
    valid = ~np.isnan(dG) & ~np.isnan(var)
    dG = dG[valid]
    var = var[valid]

    inv_var = 1.0 / var
    weighted_mean = np.sum(inv_var * dG) / np.sum(inv_var)
    weighted_var = 1.0 / np.sum(inv_var)

    return weighted_mean, weighted_var

def compute_delta_delta_G(dG_list_1, var_list_1, dG_list_2, var_list_2):
    G1, var_G1 = compute_weighted_average(dG_list_1, var_list_1)
    G2, var_G2 = compute_weighted_average(dG_list_2, var_list_2)

    delta_delta_G = G2 - G1
    var_delta_delta_G = var_G1 + var_G2
    sigma_delta_delta_G = np.sqrt(var_delta_delta_G)

    return {
        "ΔΔG": delta_delta_G,
        "var(ΔΔG)": var_delta_delta_G,
        "σ(ΔΔG)": sigma_delta_delta_G,
        "ΔG_1": G1,
        "var(ΔG_1)": var_G1,
        "ΔG_2": G2,
        "var(ΔG_2)": var_G2
    }


proj_dir = '/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45'
system_name = 'complex'
list_clones_1t = [0,1,2,3,4]
list_clones_no1t = [5,6,7,8,9]
#system_name = 'water'
#list_clones_1t = [0,1,3,4]
#list_clones_no1t = [2,5]
#swaps_per_ns, lambdas, eq_dist, k = scrape_mdp(mdp_file = f'{proj_dir}/{sys}/CLONE{c}/frame1.mdp')

free_energies_1t, swaps_per_ns_1t = load_clones(list_clones_1t)
free_energies_no1t, swaps_per_ns_no1t = load_clones(list_clones_no1t)

#swaps_per_ns = swaps_per_ns_1t  # safe to use now
#n_steps = len(free_energies_1t[0])
#time_ns = np.arange(n_steps) / swaps_per_ns

# Align all trajectories to the shortest length
#min_length = min([len(f) for f in free_energies_1t + free_energies_no1t])

#free_energies_1t = [f[:min_length] for f in free_energies_1t]
#free_energies_no1t = [f[:min_length] for f in free_energies_no1t]
#time_ns = time_ns[:min_length]

dG_list_1, var_list_1 = cal_per_clone_avg_dG(proj_dir, system_name, list_clones_1t)
weighted_avg = compute_weighted_mean_dG(dG_list_1, var_list_1)
print(f"Weighted average ΔG = {weighted_avg:.4f}")

dG_list_2, var_list_2 = cal_per_clone_avg_dG(proj_dir, system_name, list_clones_no1t)
weighted_avg = compute_weighted_mean_dG(dG_list_2, var_list_2)
print(f"Weighted average ΔG = {weighted_avg:.4f}")

result = compute_delta_delta_G(dG_list_1, var_list_1, dG_list_2, var_list_2)

# === Print Result ===
print(f"ΔG_1 = {result['ΔG_1']:.4f} ± {np.sqrt(result['var(ΔG_1)']):.4f}")
print(f"ΔG_2 = {result['ΔG_2']:.4f} ± {np.sqrt(result['var(ΔG_2)']):.4f}")
print(f"ΔΔG = {result['ΔΔG']:.4f} ± {result['σ(ΔΔG)']:.4f} (variance = {result['var(ΔΔG)']:.4f})")


# Given values
delta_delta_G = 0.9641
delta_delta_G_err = 0.6979  # std = sqrt(0.4870)

# Plot
plt.figure(figsize=(2.5, 4))

# Horizontal reference line at ΔΔG = 0
plt.axhline(0, color='magenta', linewidth=2)

# Single ΔΔG point at x=0 with error bar
plt.errorbar(x=0, y=delta_delta_G, yerr=delta_delta_G_err,
             fmt='o', color='brown', capsize=5, markersize=10)

# Axis settings
plt.xlim(-0.5, 0.5)
plt.xticks([0], ['avg(ΔΔG)'], rotation=45)
plt.ylabel("ΔΔG (kcal/mol)")
plt.title("ejm_31 → ejm_45")

# Clean layout
plt.tight_layout()
plt.box(True)
plt.grid(False)

# Show or save
plt.show()
# plt.savefig("delta_delta_G_single.pdf", dpi=600, bbox_inches='tight')

# Show or save
plt.savefig("/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45/EE_analysis/plots/delta_delta_G_plot.pdf", format='pdf', dpi=600, bbox_inches='tight')



