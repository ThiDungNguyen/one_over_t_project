import os, glob, re, sys, subprocess
import numpy as np
from matplotlib import pyplot as plt

# Command to load GROMACS 2020.3 module
module_load_command = 'module load gromacs/2020.3'
subprocess.run(module_load_command, shell=True)

proj_path = '/home/tun50867/work/git/one_over_t_project/systems/edge_ejm_31_ejm_45'
script_path = '/home/server/git'

def get_new_index_files(proj,  group = 'none', gro = 'prod'):
    '''generate a new index file containing a group of all atoms of natural and non-natural residues'''
    
    if group != 'none':
       os.system(f"echo '{group}\nq\n' | gmx make_ndx -f {proj_path}/{proj}/{gro}.gro -o {proj_path}/{proj}/{proj}.ndx")
    
def write_mdp_file(proj):
    f = open(f'{proj_path}/nonwater.mdp','w')
    f.write('''; Run parameters
integrator             = steep
nsteps                 = 100000
emtol                  = 10
emstep                 = 0.1''')
    f.close()

def get_new_tpr(proj,sys, top = 'topol', gro = 'prod'):
    'Build a custom *.tpr for the subset of atoms  in the *.xtc trajectories'
 
    # edit nonwater top file
    if not os.path.exists(f'{proj_path}/p{proj}_tops'):
        os.system(f'mkdir {proj_path}/p{proj}_tops')
    os.system(f'ln -s /home/tun50867/work/wang_lab_colab/CODEs/amber14sb.ff  {proj_path}/p{proj}_tops')
    f = open(f'{proj_path}/GMX/CLONE0/{top}.top','r')
    f_lines = f.readlines()
    f.close()
    f_final = open(f'{proj_path}/{proj}.top','w')
    for i in range(len(f_lines)):
        if 'SOL' in f_lines[i]:
            f_final.write('')
        elif 'CL' in f_lines[i] and len(f_lines[i].split()) == 2:
            f_final.write(f_lines[i])
        elif 'NA' in f_lines[i] and len(f_lines[i].split()) == 2:
            f_final.write(f_lines[i])
        else:
            f_final.write(f_lines[i])
    f_final.close()

    # generate xtc gro file
    if not os.path.exists(f'{proj_path}/p{proj}_xtc_gros'):
        os.system(f'mkdir {proj_path}/p{proj}_xtc_gros')
    os.system(f'gmx editconf -f {proj_path}/GMX/CLONE0/{gro}.gro -o {proj_path}/{proj}_xtc_gros/{sys}_1.pdb')
    os.system(f'grep -v -e "SOL" {proj_path}/{proj}_xtc_gros/{sys}_1.pdb > {proj_path}/{proj}_xtc_gros/{sys}.pdb')
    os.system(f'gmx editconf -f {proj_path}/{proj}_xtc_gros/{sys}.pdb -o {proj_path}/{proj}_xtc_gros/{sys}.gro')
    os.system(f'rm {proj_path}/{proj}_xtc_gros/{sys}_1.pdb')
    os.system(f'rm {proj_path}/{proj}_xtc_gros/{sys}.pdb')
    
    # gmx grompp to make a fake *.tpr
    if not os.path.exists(f'{proj_path}/{proj}_tprs'):
        os.system(f'mkdir {proj_path}/{proj}_tprs')
    write_mdp_file(proj)
    os.system(f'gmx grompp -c {proj_path}/{proj}_xtc..gro -f {proj_path}/nonwater.mdp -p {proj_path}/{proj}.top -o {proj_path}/{proj}.tpr')

def make_protein_only_gro(proj,  group, gro = 'prod'):
    if not os.path.exists(f'{proj_path}/{proj}_protein_gros'):
        os.system(f'mkdir {proj_path}/{proj}_protein_gros')
    os.system(f'echo "{group} \n" | gmx trjconv -s {proj_path}/{proj}.tpr -f {proj_path}/GMX/CLONE0/{gro}.gro -n {proj_path}/{proj}.ndx -o {proj_path}/{proj}_protein.gro')

def concatenate_pbc_correction(proj, ci, cf,date, top = 'topol'):
    '''concatenate and pbc correction'''

    #for r in range(ri,rf): # for nRUNs
    #make_protein_only_gro(proj, r)
    group = 2 
    get_new_index_files(proj, group = 2, gro='nvt')
    #get_new_tpr(proj, sys, top = 'topol', gro = 'prod')
    #make_protein_only_gro(proj, sys, group = group)
    if not os.path.exists(f'{proj_path}/{proj}/{proj}-TRAJECTORIES-{date}'):
        os.system(f'mkdir {proj_path}/{proj}/{proj}-TRAJECTORIES-{date}')
    for c in range(ci,cf): #for nCLONEs
        print(f'for clone{c}') 
        run_clone_path = f'{proj_path}/{proj}/CLONE{c}'
        print('run_clone_path',run_clone_path)
        if os.path.exists(f'{run_clone_path}/frame0.gro'):
            A = glob.glob(f'{run_clone_path}/frame*.gro')
            print('A',len(A))
            traj_length = len(A)*1/1000 #convert ns to us
            print(f'{run_clone_path}')
            if len(A) > 1: # 100, jusr take traj longer than 0.5 us 
                #change name of traj file to make sure they are in order : 01 - 02 -03 ... - 10 -11 ...
                for k in range(len(A)):
                    os.system(f'cp {run_clone_path}/frame{k}.xtc {run_clone_path}/traj_0{k}.xtc')
                os.system(f'gmx trjcat -f {run_clone_path}/traj_0*.xtc -o {proj_path}/{proj}/{proj}-TRAJECTORIES-{date}/C{c}_traj.xtc')
                os.system(f"echo '{group}\n{group}\n' | gmx trjconv -f {proj_path}/{proj}/{proj}-TRAJECTORIES-{date}/C{c}_traj.xtc -s {proj_path}/{proj}/pre_EE.tpr -n {proj_path}/{proj}/{proj}.ndx -pbc mol -center -o {proj_path}/{proj}/{proj}-TRAJECTORIES-{date}/CLONE{c}_{traj_length}us.xtc")
                os.system(f'rm {proj_path}/{proj}/{proj}-TRAJECTORIES-{date}/C{c}_traj.xtc')

for s in ['water']:
    #get_new_index_files(proj ='12435', sys = s, a =13)
    #get_new_tpr(proj = '12435',sys = s, top = 'topol', gro = 'prod')
    #make_protein_only_gro(proj = '12435', sys = s, group =24)
    concatenate_pbc_correction(proj='water', ci=40,cf=50, date = '20Jun2025')
    #plot_traj_length(proj = '12436', ri =9 , rf=10, ci = 0, cf = 1000 ,date ='11Dec2023', temp = 300)
    #write_traj_length_forall_clones(proj = '16977', run =r, date = '2Jan2024', ns_per_frame = 5)
    #write_traj_length_for100th_clones(proj = '16977', run = r, date = '2Jan2024', ns_per_frame =5)
    #plot_length_dist(proj = '16977', run = 13, date = '2Jan2024', lable = '_100thclones')
    #plot_length_dist(proj, run,date)
