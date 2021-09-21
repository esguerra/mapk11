#!/bin/bash
#SBATCH -J mapk11
#SBATCH --output=mapk11_%j.out
#SBATCH --error=mapk11_%j.err
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=6
#SBATCH --time=30:00:00

module purge
source /shared/work/mesguerra/miniconda3.8/bin/activate
conda activate mapk11
#export SRUN=1  # this is to avoid having to set usesrun: true in input.yaml
hostname
python mapk11_wf.py

#python -c "import pele_platform; print('Using PELEPlatform, version', pele_platform.__version__)"





