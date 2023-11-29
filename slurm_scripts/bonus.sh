#!/bin/bash
#SBATCH --job-name=mpi_openmp # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

cd /nfsmnt/120090222/CSC4005/project3
mpirun -np 3 ./bonus 1000 1000 3 >> ./test_data/bonus.txt
# mpirun -np 20 ./bonus 1000 1000 20
# mpirun -np 40 ./bonus 1000 1000 40