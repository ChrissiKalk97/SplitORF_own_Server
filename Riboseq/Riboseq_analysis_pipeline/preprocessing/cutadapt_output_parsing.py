import sys
import re
import csv
import pandas as pd


def parse_cutadapt_block(block):
    # Extract filename
    filename_match = re.search(r'-o (\S+)|--output (\S+)', block)
    filename = filename_match.group(1) or filename_match.group(
        2) if filename_match else "unknown"
    filename = filename.split('/')[-1]
    filename = filename.split('.')[0]

    # Extract total read pairs
    total_reads = int(re.search(
        r'Total read pairs processed:\s+([\d,]+)', block).group(1).replace(',', ''))

    # Reads after filtering
    reads_written = int(re.search(
        r'Pairs written \(passing filters\):\s+([\d,]+)', block).group(1).replace(',', ''))

    # Adapter hits
    r1_adapter = re.search(
        r'Read 1 with adapter:\s+([\d,]+) \(([\d\.]+)%\)', block)
    r2_adapter = re.search(
        r'Read 2 with adapter:\s+([\d,]+) \(([\d\.]+)%\)', block)

    # Adapter base composition (Read 1)
    base_stats_match = re.search(
        r'Bases preceding removed adapters:\n\s+A: ([\d.]+)%\n\s+C: ([\d.]+)%\n\s+G: ([\d.]+)%\n\s+T: ([\d.]+)%', block)
    if base_stats_match:
        base_A, base_C, base_G, base_T = map(float, base_stats_match.groups())
    else:
        base_A = base_C = base_G = base_T = None

    return {
        "filename": filename,
        "total_read_pairs": total_reads,
        "read_pairs_written": reads_written,
        "read1_adapter_count": int(r1_adapter.group(1).replace(',', '')),
        "read1_adapter_percent": float(r1_adapter.group(2)),
        "read2_adapter_count": int(r2_adapter.group(1).replace(',', '')),
        "read2_adapter_percent": float(r2_adapter.group(2)),
        "r1_base_A": base_A,
        "r1_base_C": base_C,
        "r1_base_G": base_G,
        "r1_base_T": base_T
    }


def main():
    cutadapt_out_path = sys.argv[1]
    cutadapt_summary_csv = sys.argv[2]
    # Read the whole file
    with open(cutadapt_out_path) as f:
        content = f.read()

    # Split by "This is cutadapt" to separate samples
    blocks = re.split(r'This is cutadapt .*?\n', content)[1:]  # drop preamble

    # Parse each block
    parsed_data = [parse_cutadapt_block(block) for block in blocks]

    # Write to CSV
    with open(cutadapt_summary_csv, "w", newline="") as csvfile:
        fieldnames = list(parsed_data[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(parsed_data)

    cutadapt_out = pd.read_csv(
        cutadapt_summary_csv)
    print(cutadapt_out, flush=True)


if __name__ == '__main__':
    main()
