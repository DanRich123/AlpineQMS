#!/bin/bash
#SBATCH --job-name=qms_3d
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
# #SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=normal_%j.out
#SBATCH --error=normal_%j.err

# Print job info
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load modules
module purge
module load anaconda
module load intel
module load mkl

# Activate environment
conda activate sandbox

# Run simulation
python master.py
python plot.py