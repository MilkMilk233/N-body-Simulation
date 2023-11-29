#!/bin/bash

#SBATCH --job-name=cuda  ## Job name
#SBATCH --gres=gpu:1                ## Number of GPUs required for job execution.
#SBATCH --partition=Project           ## the partitions to run in (Debug or Project)
#SBATCH --ntasks=1                  ## number of tasks (analyses) to run
#SBATCH --gpus-per-task=1           ## number of gpus per task
#SBATCH --time=0-00:02:00           ## time for analysis (day-hour:min:sec)

## Compile the cuda script using the nvcc compiler
## You can compile your codes out of the script and simply srun the executable file.
cd /nfsmnt/120090222/CSC4005/project3
## Run the script
./cuda 1000 1000 >> ./test_data/cuda.txt
# ./cuda 10000 100
# ./cuda 10000 100
