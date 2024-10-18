import os
import re
import pandas as pd
import argparse

# Function to extract values from a log file
def extract_data(log_file):
    ti_endpoint = None
    ti_error = None
    mbar_endpoint = None
    mbar_error = None
    bar_endpoint = None
    bar_error = None

    # Read the log file
    with open(log_file, 'r') as file:
        content = file.read()

        # Extract TI endpoint difference
        ti_endpoint_match = re.search(r'Endpoint differences \(TI\)\s+([\d.-]+)', content)
        if ti_endpoint_match:
            ti_endpoint = float(ti_endpoint_match.group(1))

        # Extract TI error endpoint difference
        ti_error_match = re.search(r'TI error endpoint difference\s+([\d.]+)', content)
        if ti_error_match:
            ti_error = float(ti_error_match.group(1))

        # Extract MBAR endpoint difference
        mbar_endpoint_match = re.search(r'Endpoint differences \(MBAR\)\s+([\d.-]+)', content)
        if mbar_endpoint_match:
            mbar_endpoint = float(mbar_endpoint_match.group(1))

        # Extract MBAR error endpoint difference
        mbar_error_match = re.search(r'MBAR error endpoint difference\s+([\d.]+)', content)
        if mbar_error_match:
            mbar_error = float(mbar_error_match.group(1))

        # Extract BAR endpoint difference
        bar_endpoint_match = re.search(r'Endpoint differences \(BAR\)\s+([\d.-]+)', content)
        if bar_endpoint_match:
            bar_endpoint = float(bar_endpoint_match.group(1))

        # Extract BAR error endpoint difference
        bar_error_match = re.search(r'(?<!M)BAR error endpoint difference\s+([\d.]+|nan)', content)
        if bar_error_match:
            bar_error_value = bar_error_match.group(1)
            bar_error = float(bar_error_value) if bar_error_value != 'nan' else None

    return ti_endpoint, ti_error, mbar_endpoint, mbar_error, bar_endpoint, bar_error

# Main function to process log files and save to CSV
def process_logs_to_csv(log_dir, output_csv):
    data = []

    # Loop through all log files in the specified directory
    for filename in os.listdir(log_dir):
        if filename.endswith('.log'):
            log_file_path = os.path.join(log_dir, filename)
            ti_endpoint, ti_error, mbar_endpoint, mbar_error, bar_endpoint, bar_error = extract_data(log_file_path)
            data.append({
                'file/replicate': filename,
                'TI endpoint difference': ti_endpoint,
                'TI error': ti_error,
                'MBAR endpoint difference': mbar_endpoint,
                'MBAR error': mbar_error,
                'BAR endpoint difference': bar_endpoint,
                'BAR error': bar_error
            })

    # Create a DataFrame and save to CSV
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)

# Main entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract simulation log data and save to CSV.')
    parser.add_argument('-i', '--log_directory', type=str, help='Directory containing log files')
    parser.add_argument('-o','--output_file', type=str, help='Output CSV file name')

    args = parser.parse_args()
    process_logs_to_csv(args.log_directory, args.output_file)
    print(f'Data saved to {args.output_file}')

