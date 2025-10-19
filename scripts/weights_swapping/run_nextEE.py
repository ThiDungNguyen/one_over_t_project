import numpy as np

def parse_mdp_file(file_path):
    """Parses the MDP file for `nsteps` and `dt`."""
    nsteps, dt = None, None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('nsteps'):
                nsteps = int(line.split()[2])
            elif line.startswith('dt'):
                dt = float(line.split()[2])
    if nsteps is None or dt is None:
        raise ValueError("nsteps or dt not found in MDP file.")
    return nsteps, dt

def parse_log_file(log_file, hist_counts = True):
    """Parses the log file to extract WL increment, weights, and histogram counts. For absolute binding free energy calculation"""
    with open(log_file, 'r') as f:
        lines = f.readlines()

    # Find the last occurrence of "MC-lambda information"
    start_indices = [i for i, line in enumerate(lines) if "MC-lambda information" in line]
    if not start_indices:
        raise ValueError("MC-lambda information not found in log file.")

    start_idx = start_indices[-1]

    # Collect lines from the last block
    chunk = []
    while lines[start_idx].strip():
        chunk.append(lines[start_idx])
        start_idx += 1

    # Extract information
    chunk.pop(0)  # Remove header
    wl_increment_in_kT = float(chunk.pop(0).split()[-1])
    print('wl_increment_in_kT',wl_increment_in_kT)
    chunk.pop(0)  # Remove table header

    # Parse WL weights and histogram counts
    wang_landau_weights = np.array([float(line.split()[6]) for line in chunk])

    # Find the current lambda state
    lambda_state = next((int(line.split()[0]) - 1 for line in chunk if '<<' in line), 0)

    if hist_counts:
        histogram_counts = np.array([int(line.split()[5]) for line in chunk])
        return wl_increment_in_kT, wang_landau_weights, histogram_counts, lambda_state
    else:
        return wl_increment_in_kT, wang_landau_weights, lambda_state

def modify_mdp(prev_log, prev_mdp, new_mdp, gen, hist_counts = True):
    """
    Modifies the MDP file with new start time, WL weights, histogram counts, and lambda state.
    """
    # Parse MDP for nsteps and dt
    nsteps_per_WU, dt = parse_mdp_file(prev_mdp)
    starttime = nsteps_per_WU * gen
    time_init = starttime * dt

    print(f"Starting time: {starttime}, Time Init: {time_init}, gen: {gen}")

    # Parse log for WL data
    if hist_counts:
        wl_delta, wl_weights, histogram_counts, lambda_state = parse_log_file(prev_logi, hist_counts = True)
        # Convert lists to string format for MDP
        new_weights = ' '.join(f"{w:.5f}" for w in wl_weights)
        new_histogram = ' '.join(f"{int(h)}" for h in histogram_counts)

        # Modify lines in memory
        keys_to_update = {
            "init-step":starttime,
            "tinit":time_init,
            "gen_vel":';gen_vel',
            "gen-temp":';gen-temp',
            "gen-seed":';gen-seed',
            "continuation":'yes',
            "init_lambda_state": lambda_state,
            "init-wl-delta": wl_delta,
            "init-lambda-weights": new_weights,
            "; init-lambda-weights": new_weights,
            "init-wl-histogram-counts": new_histogram,
            "; init-wl-histogram-counts": new_histogram
        }

    else:
        wl_delta, wl_weights, lambda_state = parse_log_file(prev_log,hist_counts = False)
        # Convert lists to string format for MDP
        new_weights = ' '.join(f"{w:.5f}" for w in wl_weights)

        # Modify lines in memory
        keys_to_update = {
            "init-step":starttime,
            "tinit":time_init,
            "gen_vel":';gen_vel',
            "gen-temp":';gen-temp',
            "gen-seed":';gen-seed',
            "continuation":'yes',
            "init-lambda-state": lambda_state,
            "init-wl-delta": wl_delta,
            "init-lambda-weights": new_weights,
            "; init-lambda-weights": new_weights,
        }

    # Read original MDP content
    with open(prev_mdp, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        for key, value in keys_to_update.items():
            if line.startswith(key):
                line = f"{key} = {value}\n"
                if line.startswith(';'):
                    line = line[1:].lstrip()
        new_lines.append(line)

    # Write updated content to new MDP file
    with open(new_mdp, "w") as f:
        f.writelines(new_lines)

    print(f"New MDP file saved as: {new_mdp}")

def perform_next_gen_ee_for1overt_ejm():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -l walltime=4:00:00
#PBS -o {out_dir}/{sys_name}.out
#PBS -e {out_dir}/{sys_name}.err

cd {out_dir}
export OMP_NUM_THREADS=1

module load gcc/12.2.0
export PATH=/home/tun50867/gromacs_based2025_3rd/gromacs/build/bin:$PATH
export GMXBIN=/home/tun50867/gromacs_based2025_3rd/gromacs/build/bin
export GMX_LIB=/home/tun50867/gromacs_based2025_3rd/gromacs/build/lib
export LD_LIBRARY_PATH=$GMX_LIB:$LD_LIBRARY_PATH
export GMX_EXE=$GMXBIN/gmx

#for first ee run
#gmx grompp -f frame0.mdp -c ../nvt.gro -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx

gmx grompp -f frame{gen}.mdp -c frame{gen-1}.gro -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx 
gmx mdrun -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -v -ntmpi 24 -ntomp 1

"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)



def wait_for_files(g, clones, out_dir, poll_interval=60):
    file_paths = [f"{out_dir}/frame{g}.gro"]
    print(f'waiting for a {file_paths}')
    while not all(os.path.exists(file_path) for file_path in file_paths):
        time.sleep(poll_interval)
    print(f'DONe waiting for a {file_paths}')

def perform_next_gen_ee_for1overt_paracetamol():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q normal
#PBS -l nodes=3:ppn=28
#PBS -l walltime=12:00:00
#PBS -o {out_dir}/{sys_name}.out
#PBS -e {out_dir}/{sys_name}.err

cd {out_dir}

module load gromacs/2020.3

#for first ee run
#gmx grompp -f ../Paracetamol_coul_vdW/opt_ee_optimized.mdp -c ../npt.gro -p ../system.top -o frame{gen}.tpr -maxwarn 1 

gmx grompp -f frame{gen}.mdp -c frame{gen-1}.gro -p ../system.top -o frame{gen}.tpr -maxwarn 1
gmx mdrun -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -v -ntmpi 24 -ntomp 1
#-ntmpi 24 -ntomp 1

"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

def perform_next_gen_ee_for1overt_paracetamol_Clone10():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q normal
#PBS -l nodes=3:ppn=16
#PBS -l walltime=24:00:00
#PBS -o {out_dir}/{sys_name}.out
#PBS -e {out_dir}/{sys_name}.err

cd {out_dir}
export OMP_NUM_THREADS=1

module load gcc/12.2.0

export PATH=/home/tun50867/gromacs2025_with_mpifixed/gromacs-2025-rc/build/bin:$PATH
export GMXBIN=/home/tun50867/gromacs2025_with_mpifixed/gromacs-2025-rc/build/bin
export GMX_LIB=/home/tun50867/gromacs2025_with_mpifixed/gromacs-2025-rc/build/lib
export LD_LIBRARY_PATH=$GMX_LIB:$LD_LIBRARY_PATH
export GMX_EXE=$GMXBIN/gmx

#for first ee run
gmx grompp -f ../Paracetamol_coul_vdW/opt_ee_optimized.mdp -c ../npt.gro -p ../system.top -o frame{gen}.tpr -maxwarn 1

#gmx grompp -f frame{gen}.mdp -c frame{gen-1}.gro -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx
gmx mdrun -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -v 
# -ntmpi 12 -ntomp 1

"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)


### main ###

import sys, os, time, subprocess
sys.path.append('/home/tun50867/work/git/codes')
script_path = '/home/tun50867/work/git'
from build_sys_lambopt_eeanalysis_forall import run_ee_onowls

for gen in range(14,15):
    gen +=1 
    for clone in [9]:

        out_dir = f'/home/tun50867/work/git/one_over_t_project/systems/Paracetamol/CLONE{clone}'

        wait_for_files(gen-1, clone, out_dir)

        sys_name = f'paracetamol{clone}g{gen}'
        print('run run_next_ee')
        print('modify mdp file -------------------')

        modify_mdp(prev_log = f'{out_dir}/frame{gen-1}.log', prev_mdp = f'{out_dir}/frame{gen-1}.mdp', new_mdp = f'{out_dir}/frame{gen}.mdp',gen = gen,hist_counts = False)
        #print('start perform next ee gen ------------')
        
        perform_next_gen_ee_for1overt_paracetamol()
        #perform_next_gen_ee_for1overt_paracetamol_Clone10()
