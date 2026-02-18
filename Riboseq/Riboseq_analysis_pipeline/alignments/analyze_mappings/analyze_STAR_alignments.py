import sys
import os
import csv

# Function to parse a single Log.final.out file
def parse_star_log(file_path):
    metrics = {}
    with open(file_path, 'r') as file:
        for line in file:
            if 'Uniquely mapped reads %' in line:
                metrics['Uniquely Mapped (%)'] = line.split('|')[-1].strip()
            elif '% of reads mapped to multiple loci' in line:
                metrics['Multimapped (%)'] = line.split('|')[-1].strip()
            elif '% of reads mapped to too many loci' in line:
                metrics['Too Many Loci (%)'] = line.split('|')[-1].strip()
            elif '% of reads unmapped: too many mismatches' in line:
                metrics['Too Many Mismatches (%)'] = line.split('|')[-1].strip()
            elif '% of reads unmapped: too short' in line:
                metrics['Too Short (%)'] = line.split('|')[-1].strip()
            elif '% of reads unmapped: other' in line:
                metrics['Unmapped: Other (%)'] = line.split('|')[-1].strip()
    return metrics

# Directory containing Log.final.out files
log_dir = sys.argv[1]
output_csv = sys.argv[1] + '/' + sys.argv[2]

# Parse all Log.final.out files in the directory
summary_data = []
for log_file in os.listdir(log_dir):
    if log_file.endswith('Log.final.out'):
        file_path = os.path.join(log_dir, log_file)
        sample_name = os.path.basename(log_file).replace('Log.final.out', '').strip()
        metrics = parse_star_log(file_path)
        metrics['Sample'] = sample_name
        summary_data.append(metrics)

# Write the results to a CSV file
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['Sample', 'Uniquely Mapped (%)', 'Multimapped (%)', 
                  'Too Many Loci (%)', 'Too Many Mismatches (%)', 
                  'Too Short (%)', 'Unmapped: Other (%)']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in summary_data:
        writer.writerow(row)

print(f"Summary written to {output_csv}")