import os
import statistics
from Bio import SeqIO

raw_dir = "/scratch/sharma.d2/Sanger_project/data/raw_ab1/"
trimmed_file = "/scratch/sharma.d2/Sanger_project/output/trimmed/all_trimmed.fasta"
output_file = "/scratch/sharma.d2/Sanger_project/output/final_qc_report.csv"

# Load trimmed sequences
trimmed = {}
for record in SeqIO.parse(trimmed_file, "fasta"):
    trimmed[record.id] = len(record.seq)

# Build report
with open(output_file, "w") as out:
    out.write("Sample,Original_bp,Trimmed_bp,Bases_removed,MeanQ,MedianQ,PctQ20,PctQ30\n")

    for filename in sorted(os.listdir(raw_dir)):
        if filename.endswith(".ab1"):
            filepath = os.path.join(raw_dir, filename)
            record = SeqIO.read(filepath, "abi")
            quals = record.letter_annotations["phred_quality"]

            sample_id = filename.replace(".ab1", "")
            original = len(record.seq)
            trimmed_len = trimmed.get(sample_id, 0)
            removed = original - trimmed_len
            mean_q = round(statistics.mean(quals), 1)
            median_q = round(statistics.median(quals), 1)
            pct_q20 = round(100 * sum(1 for q in quals if q >= 20) / len(quals), 1)
            pct_q30 = round(100 * sum(1 for q in quals if q >= 30) / len(quals), 1)

            print(f"{sample_id}: {original}bp → {trimmed_len}bp | MeanQ:{mean_q} | Q20:{pct_q20}% | Q30:{pct_q30}%")
            out.write(f"{sample_id},{original},{trimmed_len},{removed},{mean_q},{median_q},{pct_q20},{pct_q30}\n")

print(f"\nReport saved to {output_file}")
