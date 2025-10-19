def perform_next_gen_ee_for1overt_A16F_CLONE0():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q big
#PBS -l nodes=1:ppn=32
#PBS -l walltime=96:00:00

cd {out_dir}
export OMP_NUM_THREADS=1

module load gcc/12.2.0
export PATH=/home/tun50867/gromacs_based2025_3rd/gromacs/build/bin:$PATH
export GMXBIN=/home/tun50867/gromacs_based2025_3rd/gromacs/build/bin
export GMX_LIB=/home/tun50867/gromacs_based2025_3rd/gromacs/build/lib
export LD_LIBRARY_PATH=$GMX_LIB:$LD_LIBRARY_PATH
export GMX_EXE=$GMXBIN/gmx

gmx grompp -f frame{gen}.mdp -t frame{gen-1}.cpt -c frame{gen-1}.gro -r ../A8_A16F_NVT_equil.gro -p ../A8_A16F.top -o frame{gen}.tpr -maxwarn 1 > grompp{gen}.log

gmx mdrun -v -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -ntmpi 32 -ntomp 1 -maxh 95.9
"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

def perform_next_gen_ee_for1overt_A16F():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -l walltime=03:00:00
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

# first time run ee:
#gmx grompp -f frame{gen}.mdp -c ../A8_A16F_NVT_equil.gro -r ../A8_A16F_NVT_equil.gro -p ../A8_A16F.top -o frame{gen}.tpr -maxwarn 1 > grompp{gen}.log

gmx grompp -f frame{gen}.mdp -t frame{gen-1}.cpt -c frame{gen-1}.gro -r ../A8_A16F_NVT_equil.gro -p ../A8_A16F.top -o frame{gen}.tpr -maxwarn 1 > grompp{gen}.log

gmx mdrun -v -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -ntmpi 16 -ntomp 1 -maxh 2.9
"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

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

def perform_next_gen_ee_for1overt_ejm_CLONE10_20(sys_name, out_dir, gen):
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q medium
#PBS -l nodes=1:ppn=16
#PBS -l walltime=5:20:00
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

gmx grompp -f frame{gen}.mdp -c frame{gen-1}.gro -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx > grompp{gen}.log
gmx mdrun -s frame{gen}.tpr -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -v -ntmpi 12 -ntomp 1

"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

def perform_next_gen_ee_for1overt_complex_CLONE1():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q big
#PBS -l nodes=1:ppn=32
#PBS -l walltime=96:00:00
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

gmx grompp -f frame{gen}.mdp -t frame{gen-1}.cpt -c frame{gen-1}.gro  -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx > grompp{gen}.log

gmx mdrun -v -s frame{gen}.tpr -c frame{gen}.gro -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -ntmpi 32 -ntomp 1 -maxh 95.9


"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

def perform_next_gen_ee_for1overt_water_CLONE0():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q big
#PBS -l nodes=1:ppn=32
#PBS -l walltime=96:00:00
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

gmx grompp -f frame{gen}.mdp -t frame{gen-1}.cpt -c frame{gen-1}.gro  -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx > grompp{gen}.log

gmx mdrun -v -s frame{gen}.tpr -c frame{gen}.gro -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -ntmpi 24 -ntomp 1 -maxh 95.9


"""
    fout = open(f'{out_dir}/frame{gen}_ee.qsub', 'w')
    fout.write(qsub_text)
    fout.close()
    print(f'Wrote: {out_dir}/frame{gen}_ee.qsub')
    post_ee_cmd = f"qsub {out_dir}/frame{gen}_ee.qsub"
    subprocess.call(post_ee_cmd, shell=True)

    # Wait for the first qsub command to complete before running the second one
    subprocess.call("wait", shell=True)

def perform_next_gen_ee_for1overt_water_CLONE4():
    """Perform an initial EE simulation, with the goal of sampling all lambdas."""
    qsub_text = f"""#!/bin/bash
#PBS -N {sys_name}
#PBS -q medium
#PBS -l nodes=1:ppn=16
#PBS -l walltime=120:00:00
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

gmx grompp -f frame{gen}.mdp -t frame{gen-1}.cpt -c frame{gen-1}.gro  -p ../top/openff-1.0.0.offxml/topol1.top -o frame{gen}.tpr -maxwarn 1 -n ../index.ndx > grompp{gen}.log

gmx mdrun -v -s frame{gen}.tpr -c frame{gen}.gro -o -deffnm frame{gen} -dhdl frame{gen}.dhdl -ntmpi 12 -ntomp 1 -maxh 95.9


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

def perform_next_gen_ee_for1overt_paracetamol_Clone10_20():
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

export PATH=/home/tun50867/gromacs/build/bin:$PATH
export GMXBIN=/home/tun50867/gromacs/build/bin
export GMX_LIB=/home/tun50867/gromacs/build/lib
export LD_LIBRARY_PATH=$GMX_LIB:$LD_LIBRARY_PATH
export GMX_EXE=$GMXBIN/gmx

#for first ee run
gmx grompp -f opt_ee_optimized.mdp -c ../npt.gro -p ../system.top -o frame{gen}.tpr -maxwarn 1

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

for gen in range(0,1):
    #gen +=1 
    for clone in [10]:

        #out_dir = f'/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45/complex/CLONE{clone}'
        #out_dir = f'/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45/water/CLONE{clone}'
        #out_dir = f'/home/tun50867/work/git/one_over_t_project/systems/A16F_relative/Simulation_Files/CLONE{clone}'
        out_dir = f'/home/tun50867/work/git/one_over_t_project/systems/Paracetamol/CLONE{clone}'

        #wait_for_files(gen-1, clone, out_dir)
        o = run_ee_onowls(gen,out_dir, script_path)

        #sys_name = f'A16F_c{clone}g{gen}'
        #sys_name = f'complexc{clone}g{gen}'
        #sys_name = f'waterc{clone}g{gen}'
        sys_name = f'paracetamol{clone}g{gen}'
        print('run run_next_ee')
        print('modify mdp file -------------------')

        #o.modify_mdp(prev_gen_mdlog = f'{out_dir}/frame{gen-1}.log', prev_mdp_file = f'{out_dir}/frame{gen-1}.mdp')
        #print('start perform next ee gen ------------')

        #perform_next_gen_ee_for1overt_ejm()
        #perform_next_gen_ee_for1overt_paracetamol()
        perform_next_gen_ee_for1overt_paracetamol_Clone10_20()
        #perform_next_gen_ee_for1overt_ejm_CLONE10_20(sys_name, out_dir, gen)
        #perform_next_gen_ee_for1overt_A16F_CLONE0()
        #perform_next_gen_ee_for1overt_A16F()
        #perform_next_gen_ee_for1overt_complex_CLONE1()
        #perform_next_gen_ee_for1overt_water_CLONE0()
        #perform_next_gen_ee_for1overt_water_CLONE4()
'''
# waterc4 : run on medium, ppn: 16, timelimit: 120, mpi: 12
# waterc0: run on big, ppn: 32m timelimit 96, mpi:24
'''
