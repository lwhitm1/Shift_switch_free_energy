import argparse
import os
import shutil
import glob

def modify_coordinates(input_dir, output_dir, cmd_dir, start_atom, end_atom, iterations, lambda_states):
    try:
        gro_file = next(glob.iglob(os.path.join(input_dir, '*.gro')))
        with open(gro_file, 'r') as f:
            lines = f.readlines()

        num_atoms_line = lines[1].strip()
        num_atoms = int(num_atoms_line)
        atom_lines = lines[2:-1]
        selected_atom_lines = atom_lines[start_atom - 1:end_atom]

        for i in range(iterations):
            modified_atom_lines = []
            for line in selected_atom_lines:
                residue_number = int(line[0:5].strip())
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                atom_number = int(line[15:20].strip())
                x = float(line[20:28].strip()) + 0.001 * (i + 1)
                y = float(line[28:36].strip())
                z = float(line[36:44].strip())

                # Format the modified line according to .gro file format
                modified_line = f"{residue_number:5d}{residue_name:<5}{atom_name:>5}{atom_number:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
                modified_atom_lines.append(modified_line)

            for lambda_state in range(lambda_states):
                lambda_dir = os.path.join(output_dir, f"modified_{i}", f"lambda_{lambda_state}")
                os.makedirs(lambda_dir, exist_ok=True)
                lambda_dir = os.path.abspath(lambda_dir)

                output_file = os.path.join(lambda_dir, f"modified_{i}.gro")
                with open(output_file, 'w') as f:
                    f.write(lines[0])  # Title line
                    f.write(f"{num_atoms:5d}\n")  # Number of atoms
                    f.writelines(atom_lines[:start_atom - 1])
                    f.writelines(modified_atom_lines)
                    f.writelines(atom_lines[end_atom:])
                    f.write(lines[-1])  # Box vectors line

                # Copy all files from the input directory to the new subdirectory, except the .gro file
                for file in os.listdir(input_dir):
                    if file != os.path.basename(gro_file):
                        shutil.copy(os.path.join(input_dir, file), lambda_dir)

                # Get file names for the topology and mdp files
                top_file = next(glob.iglob(os.path.join(input_dir, '*.top')))
                mdp_file = next(glob.iglob(os.path.join(input_dir, '*.mdp')))

                 # Modify the mdp file for the current lambda state
                lambda_mdp_file = os.path.join(lambda_dir, os.path.basename(mdp_file))
                with open(mdp_file, 'r') as f:
                    mdp_content = f.readlines()

                # Update the init-lambda-state line
                for index, line in enumerate(mdp_content):
                    if 'init-lambda-state' in line:
                        mdp_content[index] = f"init-lambda-state = {lambda_state}\n"
                        break

                with open(lambda_mdp_file, 'w') as f:
                    f.writelines(mdp_content)

                # Create a list of commands for the submission files
                commands = [
                    f"cd {lambda_dir}\n",
                    f"gmx grompp -f {lambda_mdp_file} -c {output_file} -p {top_file} -o modified_{i}_lambda_{lambda_state}.tpr\n",
                    f"gmx mdrun -s modified_{i}_lambda_{lambda_state}.tpr -dhdl dhdl_modified_{i}_lambda_{lambda_state}.xvg\n"
                ]

                # Save all commands to a single .sh file
                submission_filename = os.path.join(cmd_dir, f"modified_{i}_lambda_{lambda_state}_commands.sh")
                with open(submission_filename, 'w') as file:
                    file.write("#!/bin/bash\n")
                    file.write("\n".join(commands))

                # Make the .sh file executable
                os.chmod(submission_filename, 0o755)

        print(f"Modified .gro files and command scripts have been successfully written to {output_dir} and {cmd_dir}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify coordinates in a .gro file iteratively")
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="Input .gro file")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory for modified .gro files")
    parser.add_argument("-c", "--cmd_dir", type=str, required=True, help="Output directory for .sh files for gromacs calculations")
    parser.add_argument("-s", "--start_atom", type=int, required=True, help="Start atom number")
    parser.add_argument("-e", "--end_atom", type=int, required=True, help="End atom number")
    parser.add_argument("-n", "--iterations", type=int, required=True, help="Number of iterations")
    parser.add_argument("-l", "--lambda_states", type=int, required=True, help="Number of lambda states")

    args = parser.parse_args()

    args.input_dir = os.path.abspath(args.input_dir)

    if not os.path.isdir(args.input_dir):
        print(f"Error: The input directory {args.input_dir} does not exist.")
    elif args.start_atom < 1:
        print("Error: start_atom must be greater than or equal to 1.")
    elif args.iterations < 1:
        print("Error: iterations must be greater than or equal to 1.")
    elif args.lambda_states < 1:
        print("Error: lambda_states must be greater than or equal to 1.")
    else:
        os.makedirs(args.output_dir, exist_ok=True)
        os.makedirs(args.cmd_dir, exist_ok=True)
        modify_coordinates(args.input_dir, args.output_dir, args.cmd_dir, args.start_atom, args.end_atom, args.iterations, args.lambda_states)

