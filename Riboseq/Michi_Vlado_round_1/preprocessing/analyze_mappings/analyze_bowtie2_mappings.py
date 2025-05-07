import sys
import re
import csv

# Input and output file paths
log_file = sys.argv[1]
output_csv = sys.argv[2]

# Function to parse the log file


def parse_log_file(file_path):
    summary_data = []
    sample_data = {}
    with open(file_path, "r") as file:
        for line in file:
            # Match the sample name
            sample_match = re.match(r"^=+.*Sample (\S+) finished", line)
            if sample_match:
                # Save current sample data if it's complete
                if sample_data:
                    summary_data.append(sample_data)
                    sample_data = {}
                # Start a new sample entry
                sample_data["Sample"] = sample_match.group(1)

            # Extract relevant metrics
            if "aligned 0 times" in line:
                sample_data["Unaligned (%)"] = line.split(
                    '(')[1].split('%')[0].strip()
            elif "aligned exactly 1 time" in line:
                sample_data["Aligned Once (%)"] = line.split(
                    '(')[1].split('%')[0].strip()
            elif "aligned >1 times" in line:
                sample_data["Aligned Multiple (%)"] = line.split(
                    '(')[1].split('%')[0].strip()
            elif "overall alignment rate" in line:
                sample_data["Overall Alignment Rate (%)"] = line.split('%')[
                    0].strip()

        # Append the last sample if present
        if sample_data:
            summary_data.append(sample_data)

    # Remove duplicate entries (keep only the first instance of each sample)
    seen_samples = set()
    unique_data = []
    for entry in summary_data:
        if entry["Sample"] not in seen_samples:
            unique_data.append(entry)
            seen_samples.add(entry["Sample"])

    return unique_data

# Write parsed data to CSV


def write_to_csv(data, output_path):
    fieldnames = ["Sample", "Unaligned (%)", "Aligned Once (%)",
                  "Aligned Multiple (%)", "Overall Alignment Rate (%)"]
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)


# Main execution
parsed_data = parse_log_file(log_file)
write_to_csv(parsed_data, output_csv)

print(f"Summary written to {output_csv}")
