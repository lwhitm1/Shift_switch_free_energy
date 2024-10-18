#!/bin/bash
#SBATCH --job-name=anthracene
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

cd /gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/8_9/8_9_shift-force-switch/sim_data/sim_9
gmx grompp -f /gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/8_9/8_9_shift-force-switch/input_files/expanded.mdp -c /gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/8_9/8_9_shift-force-switch/input_files/nvt.gro -p /gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/8_9/8_9_shift-force-switch/input_files/anthracene.top -o 8_9_shift-force-switch.tpr -maxwarn 4 > step1.log 2>&1
gmx mdrun -deffnm 8_9_shift-force-switch -nt 12 > step2.log 2>&1
