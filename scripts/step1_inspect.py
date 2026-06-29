import os
import statistics
from Bio import SeqIO

input_dir = "/scratch/sharma.d2/Sanger_project/data/raw_ab1/"
output_file = "/scratch/sharma.d2/Sanger_project/output/logs/inspection_report.txt"

results = []

for filename in sorted(os.listdir(input_dir)):
    if filename.endswith(".ab1"):
        filepath = os.path.join(input_dir, filename)
        record = SeqIO.read(filepath, "abi")
        quals = record.letter_annotations["phred_quality"]
        seq = str(record.seq)

        result = {
            "file": filename,
            "length": len(seq),
            "mean_q": round(statistics.mean(quals), 1),
            "min_q": min(quals),
            "max_q": max(quals),
            "sequence_preview": seq[:50]
        }
        results.append(result)

        print(f"File: {filename}")
        print(f"  Length:   {len(seq)} bp")
        print(f"  Mean Q:   {statistics.mean(quals):.1f}")
        print(f"  Min Q:    {min(quals)}")
        print(f"  Max Q:    {max(quals)}")
        print(f"  Sequence: {seq[:50]}...")
        print()

# Save report
with open(output_file, "w") as f:
    f.write("File\tLength\tMeanQ\tMinQ\tMaxQ\tPreview\n")
    for r in results:
        f.write(f"{r['file']}\t{r['length']}\t{r['mean_q']}\t{r['min_q']}\t{r['max_q']}\t{r['sequence_preview']}\n")

print(f"Report saved to {output_file}")
