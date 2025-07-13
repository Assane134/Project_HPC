#!/bin/bash
#SBATCH --job-name=omp_bench         
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=32            
#SBATCH --time=00:30:00               
#SBATCH --output=outputs/bench_openmp%j.log
#SBATCH --error=errors/bench_openmp%j.err

#SBATCH --account=HE_CS2021-UM6P-ST-SCCS-KFCI3I7PL8Q-DEFAULT-CPU

module load gcc                       

gcc -fopenmp dam_break_rusanov.c -o ./executables/dam_break_seq -lm

gcc -fopenmp dam_break_rusanov_openmp.c -o ./executables/dam_break_omp -lm

echo "threads,time" > openmp_timings.csv

# 1) Sequential run
t=$(./executables/dam_break_seq 2>&1 | awk '/Total runtime/ {print $3}')
echo "1,$t" >> openmp_timings.csv

# 2) OpenMP runs
for T in 2 4 8 16 24 32; do
  export OMP_NUM_THREADS=$T
  t=$(./executables/dam_break_omp 2>&1 | awk '/Total runtime/ {print $3}') # Extract the runtime
  echo "$T,$t" >> openmp_timings.csv
done

cat openmp_timings.csv