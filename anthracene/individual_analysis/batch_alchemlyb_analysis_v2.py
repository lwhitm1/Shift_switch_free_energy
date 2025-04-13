import os
import csv
import argparse
import numpy as np
import pandas as pd

from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.preprocessing.subsampling import slicing, decorrelate_dhdl, decorrelate_u_nk
from alchemlyb.estimators import TI, MBAR, BAR
from alchemlyb.postprocessors.units import to_kT, to_kJmol, to_kcalmol

def analyze_xvg(xvg_file, temp):
    """Run alchemlyb analysis on a single .xvg file and return estimators in all units."""
    dhdl_df = extract_dHdl(xvg_file, T=temp, filter=True)
    unk_df = extract_u_nk(xvg_file, T=temp, filter=True)

    dhdl_df = slicing(dhdl_df, lower=100)
    unk_df = slicing(unk_df, lower=100)

    dhdl_df = decorrelate_dhdl(dhdl_df, drop_duplicates=True, sort=True)
    unk_df = decorrelate_u_nk(unk_df, method='all', drop_duplicates=True, sort=True)

    def wrap_converted(delta_f, d_delta_f):
        class EstimatorResult:
            def __init__(self, delta_f, d_delta_f):
                self.delta_f_ = delta_f
                self.d_delta_f_ = d_delta_f
        return EstimatorResult(delta_f, d_delta_f)

    def extract_all(estimator_class, data):
        est = estimator_class().fit(data)
        return {
            'est': est,
            'kT': wrap_converted(to_kT(est.delta_f_, T=temp), to_kT(est.d_delta_f_, T=temp)),
            'kJ/mol': wrap_converted(to_kJmol(est.delta_f_, T=temp), to_kJmol(est.d_delta_f_, T=temp)),
            'kcal/mol': wrap_converted(to_kcalmol(est.delta_f_, T=temp), to_kcalmol(est.d_delta_f_, T=temp)),
        }

    return {
        'TI': extract_all(TI, dhdl_df),
        'MBAR': extract_all(MBAR, unk_df),
        'BAR': extract_all(BAR, unk_df)
    }

def propagate_bar_error(error_matrix):
    """Compute propagated BAR error from tridiagonal elements."""
    if error_matrix is None or error_matrix.empty:
        return None

    tridiag_errors = []
    for i in range(len(error_matrix) - 1):
        err = error_matrix.iloc[i, i + 1]
        if not np.isnan(err):
            tridiag_errors.append(err)

    if not tridiag_errors:
        return None

    tridiag_errors = np.array(tridiag_errors)
    return np.sqrt(np.sum(tridiag_errors**2))

def extract_values(estimator_dict):
    """Extract endpoint ΔG and uncertainties from estimators in all units."""
    values = {}

    for method, result in estimator_dict.items():
        for unit in ['kT', 'kJ/mol', 'kcal/mol']:
            est = result[unit]
            delta = est.delta_f_.loc[0.0, 1.0]

            if method == 'BAR':
                error = propagate_bar_error(est.d_delta_f_)
            else:
                error = est.d_delta_f_.loc[0.0, 1.0]

            # Format keys like: "TI (kJ/mol)", "TI Error (kJ/mol)"
            label = f"{method} ({unit})"
            values[label] = delta
            values[f"{method} Error ({unit})"] = error

    return values

def process_all_replicates(parent_dir, xvg_filename, temp):
    all_data = []
    for subdir in sorted(os.listdir(parent_dir)):
        sim_path = os.path.join(parent_dir, subdir)
        if os.path.isdir(sim_path):
            xvg_path = os.path.join(sim_path, xvg_filename)
            if not os.path.isfile(xvg_path):
                print(f"[Warning] Missing {xvg_filename} in {subdir}, skipping.")
                continue
            print(f"Processing {subdir}...")
            try:
                results = analyze_xvg(xvg_path, temp)
                row = extract_values(results)
                row['Replicate'] = subdir
                row['File'] = xvg_filename
                all_data.append(row)
            except Exception as e:
                print(f"[Error] Failed to analyze {subdir}: {e}")
    return all_data

def write_results_to_csv(data, output_csv):
    if not data:
        print("[Error] No data to write.")
        return

    methods = ['TI', 'MBAR', 'BAR']
    units = ['kT', 'kJ/mol', 'kcal/mol']

    ordered_columns = ['Replicate', 'File']
    for method in methods:
        for unit in units:
            ordered_columns.append(f"{method} ({unit})")
            ordered_columns.append(f"{method} Error ({unit})")

    existing_keys = set(data[0].keys())
    missing_keys = [k for k in existing_keys if k not in ordered_columns]
    full_columns = ordered_columns + sorted(missing_keys)

    avg_row = {'Replicate': 'Average', 'File': ''}
    for field in full_columns[2:]:
        values = [row.get(field) for row in data if row.get(field) is not None]
        avg_row[field] = np.mean(values) if values else None
    data.append(avg_row)

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=full_columns)
        writer.writeheader()
        writer.writerows(data)

    print(f"\n✅ Results written to: {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Batch alchemlyb analysis with cleaner CSV formatting.")
    parser.add_argument("-p", "--parent_dir", required=True, help="Parent directory containing replicate folders.")
    parser.add_argument("-x", "--xvg_filename", required=True, help="Name of the .xvg file in each replicate.")
    parser.add_argument("-o", "--output_csv", default="alchemlyb_analysis.csv", help="Output CSV filename.")
    parser.add_argument("-t", "--temp", type=float, default=300.0, help="Simulation temperature in Kelvin.")
    args = parser.parse_args()

    data = process_all_replicates(args.parent_dir, args.xvg_filename, args.temp)
    write_results_to_csv(data, args.output_csv)

if __name__ == "__main__":
    main()

