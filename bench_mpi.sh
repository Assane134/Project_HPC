#!/bin/bash
#SBATCH --job-name=mpi_bench                        
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=32                
#SBATCH --output=outputs/bench_mpi%j.log 
#SBATCH --error=errors/bench_mpi%j.err

#SBATCH --account=HE_CS2021-UM6P-ST-SCCS-KFCI3I7PL8Q-DEFAULT-CPU

module load gcc
module load intel                       

mpicc dam_break_rusanov_mpi.c -o ./executables/dam_break_mpi -lm

rm -f strong_scaling_mpi.csv
echo "procs,time" > strong_scaling_mpi.csv

# Run strong scaling
for P in 1 2 4 8 16 24 32; do 
  t=$(mpirun -np $P ./executables/dam_break_mpi 2>&1 \
      | awk '/Total runtime/ {print $3; exit}') # Extract the runtime
  echo "$P,$t" >> strong_scaling_mpi.csv
done

# Run weak scaling
nc0=50000 # Initial number of cells
rm -f weak_scaling_mpi.csv
echo "procs,NX,time" > weak_scaling_mpi.csv

for P in 1 2 4 8 16 24 32; do
  NX=$((nc0*P+1))
  echo "Running with $P processes and NX=$NX"
  t=$(
    mpirun -np $P ./executables/dam_break_mpi $NX \
      2>&1 | awk '/Total runtime/ {print $3; exit}' # Extract the runtime
  )
  echo "$P,$NX,$t" >> weak_scaling_mpi.csv
done
