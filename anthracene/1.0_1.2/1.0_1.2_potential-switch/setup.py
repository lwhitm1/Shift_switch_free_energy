"""Script is to setup anthracene simulations
python setup.py -n 5 -p sim_data -i input_files/ -o md_cmds --deffnm anthracene_15l -gmx /projects/liwh2139/software/gromacs/pkgs/gromacs-2024.2/bin/GMXRC --nt 12
"""
import os
import shutil
import argparse

#argparse setup
def parse_arguments():
    parser = argparse.ArgumentParser(description="Setup simulations for anthracene")
    parser.add_argument('--num_sims', '-n', required=True, help='Number of simulations to run')
    parser.add_argument('--parent_dir', '-p', required=False, help='Path to the parent directory for simulations', default='sim_data')
    parser.add_argument('--input_dir', '-i', required=True, help='Path to directory containing GROMACS input files.')
    parser.add_argument('--output_cmd_directory', '-o', required=False, help='Path to the output directory for .sh files', default='md_cmds')
    parser.add_argument('--deffnm', required=True, help='Set default name for GROMCAS output files')
    parser.add_argument('--gmx_executable', '-gmx', required=True, help='Path to the GROMACS executable you want to use')
    parser.add_argument('--nt', required=False, help='Number of threads to use for GROMACS')
    return parser.parse_args()

args = parse_arguments()

# Setup variables
num_sims = int(args.num_sims)
parent_dir = os.path.abspath(args.parent_dir)
input_dir = os.path.abspath(args.input_dir)
output_cmd_directory = os.path.abspath(args.output_cmd_directory)
gmx_executable = os.path.abspath(args.gmx_executable)
deffnm = args.deffnm
nt = args.nt

# Create directories for simulations
if not os.path.exists(parent_dir):
    os.mkdir(parent_dir)

if not os.path.exists(output_cmd_directory):
    os.mkdir(output_cmd_directory)

# Create subdirectories & copy input files for each simulation
for i in range(num_sims):
    sim_dir = os.path.join(parent_dir, f'sim_{i}')
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    # Get paths for input files & copy them to the simulation directory
    input_files = os.listdir(input_dir)
    gro_file = os.path.join(input_dir, [file for file in os.listdir(input_dir) if file.endswith('.gro')][0])
    top_file = os.path.join(input_dir, [file for file in os.listdir(input_dir) if file.endswith('.top')][0])
    mdp_file = os.path.join(input_dir, [file for file in os.listdir(input_dir) if file.endswith('.mdp')][0]) 
    
    shutil.copy2(gro_file, sim_dir)
    shutil.copy2(top_file, sim_dir)
    shutil.copy2(mdp_file, sim_dir)


    # Generate commands for each simulation
    commands = [
        "#SBATCH --job-name=anthracene\n"
        "#SBATCH --output=anthracene.%j.out\n"
        "#SBATCH --error=anthracene.%j.err\n"
        "#SBATCH --account=ucb-general\n"
        "#SBATCH --partition=amilan\n"
        "#SBATCH --qos=long\n"
        "#SBATCH --time=48:00:00\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --constraint=ib\n"
        "#SBATCH --exclude=c3cpu-c11-u19-1,c3cpu-c11-u28-2,c3cpu-c11-u34-2,c3cpu-a9-u3-2,c3cpu-c11-u15-3,c3cpu-c11-u11-4,c3cpu-c11-u13-2\n"
        f"#SBATCH --ntasks={nt}\n"
        "\n"
        "module purge\n"
        "ml gcc\n"
        "ml openmpi/5.0.6\n"
        "ml anaconda\n"
        'export OMPI_MCA_btl="self,openib,vader,tcp"\n'
        'export OMPI_MCA_pml="ob1"\n'
        "conda activate openff-toolkit\n"
        f"source {gmx_executable}\n"
        "\n"
        f"cd {sim_dir}\n"
        f"gmx grompp -f {mdp_file} -c {gro_file} -p {top_file} -o {deffnm}.tpr -maxwarn 4 > step1.log 2>&1\n"
        f"gmx mdrun -deffnm {deffnm} -nt {nt} > step2.log 2>&1\n"
    ]

    # Write commands to output file
    output_file = os.path.join(output_cmd_directory, f"sim_{i}.sh")
    with open(output_file, 'w') as file:
        file.write('#!/bin/bash\n')
        file.write('\n'.join(commands))

    # Make the .sh file executable
    os.chmod(output_file, 0o755)

# Print a message indicating the process is complete
print("Commands saved to .sh files in the output directory.")
    
