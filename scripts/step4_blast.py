from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

fasta_file = "/scratch/sharma.d2/Sanger_project/output/trimmed/all_trimmed.fasta"
output_file = "/scratch/sharma.d2/Sanger_project/output/alignments/blast_results.txt"

records = list(SeqIO.parse(fasta_file, "fasta"))

with open(output_file, "w") as out:
    out.write("Sample\tTop_Hit\tIdentity%\tE-value\tScore\n")

    for record in records:
        print(f"BLASTing {record.id}... (this may take 1-2 mins)")

        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=str(record.seq)
        )

        blast_record = NCBIXML.read(result_handle)

        if blast_record.alignments:
            top = blast_record.alignments[0]
            hsp = top.hsps[0]
            identity = round(100 * hsp.identities / hsp.align_length, 1)
            evalue = hsp.expect
            score = hsp.score
            title = top.title[:80]

            print(f"  Hit:      {title}")
            print(f"  Identity: {identity}%")
            print(f"  E-value:  {evalue}")
            print(f"  Score:    {score}")
            print()

            out.write(f"{record.id}\t{title}\t{identity}\t{evalue}\t{score}\n")
        else:
            print(f"  No hits found")
            out.write(f"{record.id}\tNo hits\tN/A\tN/A\tN/A\n")

print(f"Results saved to {output_file}")
