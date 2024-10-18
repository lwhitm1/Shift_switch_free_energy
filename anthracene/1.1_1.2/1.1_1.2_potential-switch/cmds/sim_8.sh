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

cd /scratch/alpine/liwh2139/more_lambdas/new_cutoffs/2024.3/1.1_1.2/1.1_1.2_potential-switch/sim_data/sim_8
gmx grompp -f /scratch/alpine/liwh2139/more_lambdas/new_cutoffs/2024.3/1.1_1.2/1.1_1.2_potential-switch/input_files/expanded.mdp -c /scratch/alpine/liwh2139/more_lambdas/new_cutoffs/2024.3/1.1_1.2/1.1_1.2_potential-switch/input_files/nvt.gro -p /scratch/alpine/liwh2139/more_lambdas/new_cutoffs/2024.3/1.1_1.2/1.1_1.2_potential-switch/input_files/anthracene.top -o 1.1_1.2_potential-switch.tpr -maxwarn 4 > step1.log 2>&1
gmx mdrun -deffnm 1.1_1.2_potential-switch -nt 12 > step2.log 2>&1
