#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -N Paracetamol
#PBS -q normal 
#PBS -l nodes=3:ppn=28
#PBS
cd /home/tun50867/work/git/one_over_t_project/systems/Paracetamol

module load gromacs/2020.3

# note --using ee_test.gro as starting structure since it's compatible with lam_k=0 ensemble
gmx grompp -f opt_ee_optimized.mdp -c npt.gro -r npt.gro -p system.top -o pre_ee.tpr -v -maxwarn 2
gmx mdrun -v -s pre_ee.tpr -c pre_ee.gro -g pre_ee.log -o pre_ee.trr -x pre_ee.xtc -dhdl pre_ee.dhdl 
