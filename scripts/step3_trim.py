import os
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

input_dir = "/scratch/sharma.d2/Sanger_project/data/raw_ab1/"
output_file = "/scratch/sharma.d2/Sanger_project/output/trimmed/all_trimmed.fasta"
quality_threshold = 20
min_length = 100

results = []

for filename in sorted(os.listdir(input_dir)):
    if filename.endswith(".ab1"):
        filepath = os.path.join(input_dir, filename)
        record = SeqIO.read(filepath, "abi")
        quals = record.letter_annotations["phred_quality"]
        seq = str(record.seq)

        # Find start — first position where quality exceeds threshold
        start = 0
        for i, q in enumerate(quals):
            if q >= quality_threshold:
                start = i
                break

        # Find end — last position where quality exceeds threshold
        end = len(quals)
        for i in range(len(quals) - 1, -1, -1):
            if quals[i] >= quality_threshold:
                end = i + 1
                break

        trimmed_seq = seq[start:end]
        trimmed_len = len(trimmed_seq)

        print(f"{filename}")
        print(f"  Original:  {len(seq)} bp")
        print(f"  Trimmed:   {trimmed_len} bp")
        print(f"  5' cut:    {start} bases")
        print(f"  3' cut:    {len(seq) - end} bases")
        print()

        if trimmed_len >= min_length:
            results.append(SeqRecord(
                seq=record.seq[start:end],
                id=filename.replace(".ab1", ""),
                description="trimmed"
            ))

# Save all trimmed sequences to one fasta file
SeqIO.write(results, output_file, "fasta")
print(f"Saved {len(results)} trimmed sequences to {output_file}")
