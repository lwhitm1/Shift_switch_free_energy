import os
import pandas as pd
import argparse
import re

# Function to parse a single xvg file
def parse_xvg(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    data = []
    for line in lines:
        if not line.startswith(('@', '#')):
            data.append([float(x) for x in line.split()])
    
    return data

def main(input_directory, output_file):
    # Collect all data
    all_data = []

    # Iterate over all configuration directories
    for config_dir in os.listdir(input_directory):
        if config_dir.startswith("modified_"):
            config_number = int(config_dir.split('_')[1])
            distance = (config_number + 200) / 1000.0
            
            config_path = os.path.join(input_directory, config_dir)
            
            # Iterate through each lambda directory
            for lambda_dir in os.listdir(config_path):
                if lambda_dir.startswith("lambda_"):
                    lambda_state = int(lambda_dir.split('_')[1])
                    
                    lambda_path = os.path.join(config_path, lambda_dir)
                    
                    # Locate the xvg file
                    for xvg_file in os.listdir(lambda_path):
                        if xvg_file.endswith(".xvg"):
                            xvg_path = os.path.join(lambda_path, xvg_file)
                            
                            # Parse the xvg file
                            data = parse_xvg(xvg_path)
                            
                            # Append the data to the list with additional info
                            for row in data:
                                all_data.append([xvg_file, distance, lambda_state] + row)

    # Define the columns based on the expected data structure
    columns = ["XVG Filename", "Distance (nm)", "Lambda State", "Time (ps)", 
               "Total Energy (kJ/mol)", "dH/dλ vdw-lambda = 0.0000", "ΔH λ to 0.0000", "ΔH λ to 0.1000",
               "ΔH λ to 0.2000", "ΔH λ to 0.3000", "ΔH λ to 0.4000", "ΔH λ to 0.5000",
               "ΔH λ to 0.6000", "ΔH λ to 0.7000", "ΔH λ to 0.8000", "ΔH λ to 0.9000",
               "ΔH λ to 1.0000"]

    # Convert the list to a DataFrame
    df = pd.DataFrame(all_data, columns=columns)

    # Sort the DataFrame by Distance and Lambda State
    df = df.sort_values(by=["Distance (nm)", "Lambda State"])

    # Write the DataFrame to a CSV file
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert multiple xvg files to a single CSV file.")
    parser.add_argument("-i", "--input_directory", help="Parent directory containing the configuration subdirectories")
    parser.add_argument("-o", "--output_file", help="Name of the output CSV file")

    args = parser.parse_args()
    main(args.input_directory, args.output_file)

