import os
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

trimmed_file = "/scratch/sharma.d2/Sanger_project/output/trimmed/all_trimmed.fasta"
output_file = "/scratch/sharma.d2/Sanger_project/output/alignments/reconciliation.txt"

# Load all trimmed sequences into a dictionary
sequences = {}
for record in SeqIO.parse(trimmed_file, "fasta"):
    sequences[record.id] = record.seq

# Find forward/reverse pairs
samples = set()
for name in sequences.keys():
    sample = name.replace("_F", "").replace("_R", "")
    samples.add(sample)

with open(output_file, "w") as out:
    out.write("Sample\tF_length\tR_length\tAlignment_score\tAgreement%\n")

    for sample in sorted(samples):
        f_key = sample + "_F"
        r_key = sample + "_R"

        if f_key not in sequences or r_key not in sequences:
            print(f"Skipping {sample} — missing F or R read")
            continue

        forward = sequences[f_key]
        reverse = sequences[r_key].reverse_complement()

        print(f"Sample: {sample}")
        print(f"  Forward:            {len(forward)} bp")
        print(f"  Reverse (RC):       {len(reverse)} bp")

        # Pairwise alignment
        alignments = pairwise2.align.globalms(
            str(forward), str(reverse),
            2, -1, -2, -0.5
        )

        best = alignments[0]
        q_aligned = best[0]
        r_aligned = best[1]

        # Count matches, mismatches, gaps
        matches = 0
        mismatches = 0
        gaps = 0
        variants = []

        for i, (q, r) in enumerate(zip(q_aligned, r_aligned)):
            if q == '-' or r == '-':
                gaps += 1
            elif q == r:
                matches += 1
            else:
                mismatches += 1
                variants.append(f"pos{i+1}:{q}vs{r}")

        agreement = round(100 * matches / len(q_aligned), 1)

        print(f"  Matches:            {matches}")
        print(f"  Mismatches:         {mismatches}")
        print(f"  Gaps:               {gaps}")
        print(f"  Agreement:          {agreement}%")
        if variants:
            print(f"  Variants:           {', '.join(variants[:5])}")
        print()

        out.write(f"{sample}\t{len(forward)}\t{len(reverse)}\t{round(best[2],1)}\t{agreement}\n")

print(f"Reconciliation saved to {output_file}")
