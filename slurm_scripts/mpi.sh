#!/bin/bash
#SBATCH --job-name=mpi # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

cd /nfsmnt/120090222/CSC4005/project3
mpirun -np 3 ./mpi 1000 1000 >> ./test_data/mpi.txt
# mpirun -np 20 ./mpi 1000 1000
# mpirun -np 40 ./mpi 1000 1000