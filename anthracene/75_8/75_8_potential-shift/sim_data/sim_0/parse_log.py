import re
import csv

# Define the input log file and output CSV file
log_file = '/gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/2024.3/75_8/75_8_force_switch-potential-shift/sim_data/sim_0/75_8-potential-shift.log'  # Change this to your log file name
output_csv = 'lambda_states.csv'
lambda_data = []

# Open the log file and read through it
with open(log_file, 'r') as file:
    for line in file:
        # Check for the start of the MC-lambda information block
        if "MC-lambda information" in line:
            # Read lines until we reach the lambda lines
            while True:
                line = next(file)
                # Look for the lambda line with the '<<' indicator
                if '<<' in line:
                    match_lambda = re.search(r'(\d+)\s+([\d.]+)\s*<<', line)
                    if match_lambda:
                        lambda_value = match_lambda.group(2)
                        break  # We've found the lambda value, exit the loop

            # After finding the lambda value, look for the step and time lines
            while True:
                line = next(file)
                if "Step" in line and "Time" in line:
                    # Get the next line which contains the step and time values
                    next_line = next(file)
                    match_step_time = re.search(r'(\d+)\s+([\d.]+)', next_line)
                    if match_step_time:
                        current_step = match_step_time.group(1)
                        current_time = match_step_time.group(2)
                        # Store the extracted data
                        lambda_data.append((current_step, current_time, lambda_value))
                    break  # Exit after processing the step and time

# Write the extracted data to a CSV file
with open(output_csv, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    # Write the header
    csv_writer.writerow(['Step', 'Time', 'Lambda Value'])
    # Write the data
    csv_writer.writerows(lambda_data)

print(f'Lambda states extracted to {output_csv}')

