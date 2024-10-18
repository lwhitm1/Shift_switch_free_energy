#!/bin/bash
#SBATCH --job-name=restart_anthracene
#SBATCH --output=anthracene.%j.out
#SBATCH --error=anthracene.%j.err
#SBATCH --account=ucb-general
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12

module purge
ml gcc
ml openmpi/4.1.1
source /projects/liwh2139/software/gromacs/pkgs/gromacs_2024.3/bin/GMXRC

cd /gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/75_8/75_8_force_switch-potential-shift/sim_data/sim_0

gmx mdrun -s 75_8-potential-shift.tpr -cpi 75_8-potential-shift.cpt -deffnm 75_8-potential-shift -nt 12 > step_rerun.log 2>&1
