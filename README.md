# Sanger Sequencing Analysis Pipeline

I built this pipeline as part of learning bioinformatics on a real HPC cluster. 
It takes raw Sanger sequencing files (.ab1 format) and processes them all the 
way from binary trace data to identified sequences with quality reports.

Everything was built and run on Northeastern University's Discovery cluster 
using SLURM for job scheduling and conda for environment management.


## What it does

Raw Sanger sequencing gives you .ab1 files — binary files containing the 
fluorescence signal the sequencer recorded as DNA fragments passed through 
a capillary. This pipeline reads those files, checks their quality, trims 
the noisy ends, identifies what organism the sequence came from using BLAST, 
and compares the forward and reverse reads against each other to validate 
the results.

The six steps are:

1. **Inspect** - reads every .ab1 file and pulls out the sequence length, 
   mean Phred quality score, and a preview of the first 50 bases
2. **Chromatogram** - plots the raw fluorescence trace for each sample so 
   you can visually see where the signal is clean and where it degrades
3. **Trim** - removes low quality bases from both ends using a Q20 threshold, 
   keeping only the reliable middle region
4. **BLAST** - sends each trimmed sequence to NCBI and identifies what it matches
5. **Reconcile** - reverse complements the reverse reads and aligns them against 
   the forward reads to check they agree
6. **Report** - combines everything into a single CSV with all quality metrics


## The data

I used real published Sanger sequencing data of the COI gene 
(Cytochrome Oxidase I) from *Allolobophora chlorotica*, a species of 
European earthworm. COI is the standard gene used for DNA barcoding — 
identifying species from a short DNA sequence, the same way a barcode 
identifies a product in a supermarket.

The data came from a published study deposited on Zenodo.


## What I found

All six reads matched *Allolobophora chlorotica* with E-values of 0.0 
meaning the matches are statistically certain. Samples 1 and 2 matched 
one voucher specimen (AC1) at 99.5% identity, while sample 3 matched a 
different voucher (B1) at 98.7% - suggesting it represents a slightly 
different haplotype of the same species.

The reads were very high quality - median Phred score of 61 across all 
samples, with over 92% of bases above Q30. Trimming removed fewer than 
22 bases from any read, meaning almost no data was lost.

| Sample | Identity | E-value | Mean Q |
|--------|----------|---------|--------|
| sample1_F | 99.5% | 0.0 | 52.8 |
| sample1_R | 99.1% | 0.0 | 54.9 |
| sample2_F | 99.5% | 0.0 | 55.7 |
| sample2_R | 99.4% | 0.0 | 54.5 |
| sample3_F | 98.7% | 0.0 | 54.2 |
| sample3_R | 98.9% | 0.0 | 51.6 |


## How to run it

You need Python 3.10 and BioPython:

```bash
conda env create -f environment.yml
conda activate sanger_env
```

Then run each step in order:

```bash
python scripts/step1_inspect.py
python scripts/step2_chromatogram.py
python scripts/step3_trim.py
python scripts/step4_blast.py
python scripts/step5_reconcile.py
python scripts/step6_report.py
```

If you are on an HPC cluster with SLURM, the BLAST step can be submitted 
as a batch job:

```bash
sbatch scripts/job_blast.sh
```

---

## What I learned

This project taught me how Sanger sequencing data is structured at the 
file level, how Phred quality scores work mathematically, how to submit 
and chain jobs on SLURM, and how to interpret BLAST results in a 
biological context. It was the first pipeline I built end to end on 
a real cluster with real data.

---

## Author

Divya Sharma — Northeastern University
