#!/bin/bash
#SBATCH --job-name=sanger_blast
#SBATCH --output=/scratch/sharma.d2/Sanger_project/output/logs/blast_%j.out
#SBATCH --error=/scratch/sharma.d2/Sanger_project/output/logs/blast_%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

source /shared/EL9/explorer/anaconda3/2024.06/etc/profile.d/conda.sh
conda activate /scratch/sharma.d2/envs/sanger_env

cd /scratch/sharma.d2/Sanger_project/scripts/

python step4_blast.py
