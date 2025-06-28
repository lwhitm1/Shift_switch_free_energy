#!/bin/bash
#SBATCH --job-name=anthracene
#SBATCH --output=anthracene.%j.out
#SBATCH --error=anthracene.%j.err
#SBATCH --account=ucb-general
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --constraint=ib
#SBATCH --exclude=c3cpu-c11-u19-1,c3cpu-c11-u28-2,c3cpu-c11-u34-2,c3cpu-a9-u3-2,c3cpu-c11-u15-3,c3cpu-c11-u11-4,c3cpu-c11-u13-2
#SBATCH --ntasks=12

module purge
ml gcc
ml openmpi/5.0.6
ml anaconda
export OMPI_MCA_btl="self,openib,vader,tcp"
export OMPI_MCA_pml="ob1"
conda activate openff-toolkit
source /projects/liwh2139/software/gromacs/pkgs/gromacs_2024.3/bin/GMXRC

cd /gpfs/alpine1/scratch/liwh2139/lennard-jones-miniproject/anthracene/new_cutoffs/paper_extra_Sims/0.9-1.1/potential_shift/sim_data/sim_4
gmx grompp -f /gpfs/alpine1/scratch/liwh2139/lennard-jones-miniproject/anthracene/new_cutoffs/paper_extra_Sims/0.9-1.1/potential_shift/input_files/expanded.mdp -c /gpfs/alpine1/scratch/liwh2139/lennard-jones-miniproject/anthracene/new_cutoffs/paper_extra_Sims/0.9-1.1/potential_shift/input_files/nvt.gro -p /gpfs/alpine1/scratch/liwh2139/lennard-jones-miniproject/anthracene/new_cutoffs/paper_extra_Sims/0.9-1.1/potential_shift/input_files/anthracene.top -o anthracene_psh_0.9_1.1.tpr -maxwarn 4 > step1.log 2>&1
gmx mdrun -deffnm anthracene_psh_0.9_1.1 -nt 12 > step2.log 2>&1
