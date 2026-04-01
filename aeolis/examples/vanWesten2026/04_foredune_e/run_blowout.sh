#!/bin/bash
#SBATCH --job-name=fd_04
#SBATCH --output=fd_04.o%j
#SBATCH --partition=4vcpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00

echo "Loading miniforge module..."
module load miniforge/latest || {
    echo "Failed to load miniforge"
    exit 1
}

echo "Initializing conda..."
source /opt/miniforge3/etc/profile.d/conda.sh || {
    echo "Failed to initialize conda"
    exit 1
}

ENV_PATH="/p/11211460-msc-gw-coupling/00_environments/aeolis_zandmotor_env"
echo "Activating conda environment: $ENV_PATH"
conda activate "$ENV_PATH" || {
    echo "Failed to activate conda environment"
    exit 1
}

RUNDIR="/p/11211460-msc-gw-coupling/03_simulations/fd_04"
echo "Changing to directory: $RUNDIR"
cd "$RUNDIR" || {
    echo "Failed to change directory"
    exit 1
}

echo "Running AeoLiS..."
aeolis run ./aeolis.txt || {
    echo "AeoLiS execution failed"
    exit 1
}

echo "Processing completed successfully."