#!/bin/bash
#SBATCH --account=def-stbil30-ab
#SBATCH --cpus-per-task=10
#SBATCH --mem 16G
#SBATCH -t 4:00:00              # time (D-HH:MM)
#SBATCH --job-name=Genome-Scale-Metabo
#SBATCH --output=gem%j.out
#SBATCH --error=gem%j.err

module load StdEnv/2020 python/3.11
virtualenv --no-download $SLURM_TMPDIR/ENV
source $SLURM_TMPDIR/ENV/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r requirements.txt

python GEM.py
