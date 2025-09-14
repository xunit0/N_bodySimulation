#!/bin/bash
#SBATCH --job-name=xM
#SBATCH --partition=Centaurus
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=nbody_%j.out
#SBATCH --error=nbody_%j.err

cd /users/xmolina/work/N_bodySimulation/
make

/usr/bin/time -v srun ./nbody sem 200 50000 1000


