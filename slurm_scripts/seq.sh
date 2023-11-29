#!/bin/bash
#SBATCH --job-name=seq # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # number of processes = 1 
#SBATCH --cpus-per-task=40      # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
#SBATCH --time=0-00:01:00           ## time for analysis (day-hour:min:sec)

cd /nfsmnt/120090222/CSC4005/project3
./seq 50 100 >> ./test_data/seq.txt
./seq 50 100 >> ./test_data/seq.txt
./seq 50 100 >> ./test_data/seq.txt
# ./seq 1000 1000
# ./seq 1000 1000
# ./seq 1000 1000
# ./seq 1000 1000
# ./seq 1000 1000