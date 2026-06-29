import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

input_dir = "/scratch/sharma.d2/Sanger_project/data/raw_ab1/"
output_dir = "/scratch/sharma.d2/Sanger_project/output/chromatograms/"

for filename in sorted(os.listdir(input_dir)):
    if filename.endswith(".ab1"):
        filepath = os.path.join(input_dir, filename)
        record = SeqIO.read(filepath, "abi")

        channels = {
            'G': record.annotations['abif_raw']['DATA9'],
            'A': record.annotations['abif_raw']['DATA10'],
            'T': record.annotations['abif_raw']['DATA11'],
            'C': record.annotations['abif_raw']['DATA12'],
        }

        colours = {'A': 'green', 'T': 'red', 'G': 'black', 'C': 'blue'}

        plt.figure(figsize=(15, 5))
        for base, trace in channels.items():
            plt.plot(trace, color=colours[base], label=base, alpha=0.7)

        plt.title(f"Chromatogram: {filename}")
        plt.xlabel("Scan position")
        plt.ylabel("Fluorescence intensity")
        plt.legend()
        plt.tight_layout()

        outfile = os.path.join(output_dir, filename.replace(".ab1", "_chromatogram.png"))
        plt.savefig(outfile, dpi=150)
        plt.close()
        print(f"Saved: {outfile}")

print("All chromatograms done.")
