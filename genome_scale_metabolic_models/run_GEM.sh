#!/bin/bash
#SBATCH --account=def-stbil30-ab
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --job-name=Genome-Scale-Metabo
#SBATCH --output=gem-%j.out
#SBATCH --error=gem-%j.err

# --- 1. Modules --------------------------------------------------------------
module load StdEnv/2023 python/3.12   # <- version supportée
export OMP_NUM_THREADS=1               # coupe le multithreading BLAS
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- 2. Virtualenv -----------------------------------------------------------
VENV=$SLURM_TMPDIR/env
python -m venv $VENV           # --no-download n’est plus nécessaire avec venv
source $VENV/bin/activate

# Si tu as préparé un wheelhouse local, décompresse-le ici
# tar -xzf $PROJECT/wheels.tar.gz -C $SLURM_TMPDIR/wheels

pip install --no-index --upgrade pip
pip install --no-index -r $PROJECT/genome_scale_metabolic_models/requirements.txt \
            # --find-links $SLURM_TMPDIR/wheels   # si wheelhouse perso

# --- 3. Lancement ------------------------------------------------------------
cd   $PROJECT/genome_scale_metabolic_models
python GEM.py
