import os
import csv
import re
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

from_folder = 'D:/pythonProject1/AVE'
result_folder = 'D:/pythonProject1/AVE'

if not os.path.exists(result_folder):
    os.makedirs(result_folder)

def get_all_dna_files():
    return [x for x in os.listdir(from_folder) if x.endswith('.fasta') or x.endswith('.fa') or x.endswith('.gb') or x.endswith('.gbk')]

def reverse_sequence(sequence):
    return str(Seq(sequence).reverse_complement())

def find_first_non_n_sequence(sequence):
    match = re.search(r'NN', sequence)
    if match:
        start_pos = match.end()
        while start_pos < len(sequence) and sequence[start_pos] == 'N':
            start_pos += 1
        return sequence[start_pos:]
    return sequence

def find_patterns(sequence, pattern, length, mode):
    matches = []
    sequence_len = len(sequence)
    for match in re.finditer(pattern, sequence):
        start, end = None, None
        if mode == "end":
            start = match.start() - (length - len(match.group()))
            end = match.end()
        elif mode == "start":
            start = match.start()
            end = match.start() + length
        if start is not None and end is not None and start >= 0 and end <= sequence_len:
            matches.append(sequence[start:end])
    return matches

def write_results_to_csv(result_folder, file_name, pattern, results, length):
    result_file_path = os.path.join(result_folder, f"{file_name}_{pattern}_{length}_result.csv")
    with open(result_file_path, 'w', newline='', encoding='utf8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence", "Count"])
        for result, count in results.items():
            writer.writerow([result, count])

def summarize_results(result_folder, file_name, summary):
    result_file_path = os.path.join(result_folder, f"{file_name}_result_summary.csv")
    with open(result_file_path, 'w', newline='', encoding='utf8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Pattern", "Length", "Count"])
        for pattern, length_counts in summary.items():
            for length, count in length_counts.items():
                writer.writerow([pattern, length, count])

def process_each_file(file_name, patterns, all_summaries):
    file_path = os.path.join(from_folder, file_name)
    summary = defaultdict(lambda: defaultdict(int))
    result_subfolder = os.path.join(result_folder, f"{file_name}_result")

    if not os.path.exists(result_subfolder):
        os.makedirs(result_subfolder)

    if file_name.endswith('.fasta') or file_name.endswith('.fa'):
        format = "fasta"
    elif file_name.endswith('.gb') or file_name.endswith('.gbk'):
        format = "genbank"
    else:
        print(f"Unsupported file format: {file_name}")
        return

    for record in SeqIO.parse(file_path, format):
        sequence = str(record.seq).upper()
        sequence = find_first_non_n_sequence(sequence)
        reverse_seq = reverse_sequence(sequence)

        if not sequence:
            print(f"No sequence found in {file_name}")
            continue

        print(f"Processing record: {record.id} with length {len(sequence)}")

        for pattern, mode in patterns:
            for length in range(20, 24):#12,24
                results = []
                if mode == "end":
                    regex_pattern = re.escape(pattern)
                else:
                    regex_pattern = f"(?={re.escape(pattern)})"

                # Find patterns in both forward and reverse sequences
                results.extend(find_patterns(sequence, regex_pattern, length, mode))
                results.extend(find_patterns(reverse_seq, regex_pattern, length, mode))

                results_count = Counter(results)
                filtered_results = {seq: count for seq, count in results_count.items() if count >= 5}

                if filtered_results:
                    write_results_to_csv(result_subfolder, file_name, f"{pattern}_{mode}", filtered_results, length)
                    summary[f"{pattern}_{mode}"][length] = len(filtered_results)

                # Debug output
                print(f"Pattern: {pattern}, Mode: {mode}, Length: {length}, Matches found: {len(results)}, Filtered: {len(filtered_results)}")

    summarize_results(result_subfolder, file_name, summary)
    for pattern, length_counts in summary.items():
        for length, count in length_counts.items():
            all_summaries.append([f"{file_name}_{pattern}", length, count])
    print(f"Processed file: {file_name}")


patterns = [
    #("GG", "end"),  # For GG pattern, it should be at the end
   # ("TT", "start"),  # For TT pattern, it should be at the start
    ("GG", "end"),
    # Add more patterns as needed
]

all_summaries = []

for each in get_all_dna_files():
    print(f"Start processing file: {each}")
    process_each_file(each, patterns, all_summaries)

# Write the overall summary
summary_all_path = os.path.join(result_folder, "summary_all.csv")
with open(summary_all_path, 'w', newline='', encoding='utf8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Pattern", "Length", "Count"])
    for summary in all_summaries:
        writer.writerow(summary)
