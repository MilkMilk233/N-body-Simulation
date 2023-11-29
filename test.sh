make clean
make all
sbatch ./slurm_scripts/bonus.sh
sbatch ./slurm_scripts/mpi.sh
sbatch ./slurm_scripts/cuda.sh
sbatch ./slurm_scripts/bonus.sh
sbatch ./slurm_scripts/openmp.sh
./seq 1000 1000 >> ./test_data/seq.txt