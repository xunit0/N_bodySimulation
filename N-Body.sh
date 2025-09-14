#!/bin/bash
#SBATCH --job-name=xM
#SBATCH --partition=Centaurus
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --mem=16G


cd /users/xmolina/work/N_bodySimulation/
make

/usr/bin/time -v srun ./nbody sem 200 5000000 1000


