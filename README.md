# Projet HPC: dam break problem
## Projects
### OpenMP implementation
- parallelize the loops in the code
### MPI implementation
- create many processes
- divide the data and each process will work on a part of the data
- the processes will communicate
- at the end, the results will be gathered

## Quick start
### Run OpenMP implementation
To compile the code, run:
```bash
gcc -fopenmp dam_break_rusanov_openmp.c -o ./executables/dam_break_omp -lm
```
And then, you can run the code with for example 4 threads:
```bash
export OMP_NUM_THREADS=4
./executables/dam_break_omp
```
The output file of the code is `output_openmp.dat`. With that file, you can plot the u and h curves using the python code `plot_results.py`

You can also run directly the code using the slurm job file `bench_openmp.sh`:
```bash
sbatch bench_openmp.sh
```
This way, you will launch the code with different number of threads. This will write the runtime of each number of threads in the file `openmp_timings.csv`. You can then plot the results (speedup and efficiency) with the python code `plot_openmp.py`
### Run MPI implementation
To compile the code, run:
```bash
mpicc dam_break_rusanov_mpi.c -o ./executables/dam_break_mpi -lm
```
And then, you can run the code with for example 4 processes:

```bash
mpirun -np 4 ./executables/dam_break_mpi
```
The output file of the code is `output_mpi.dat`. With that file, you can plot the u and h curves using the python code `plot_results.py`

You can also run directly the code using the slurm job file `bench_mpi.sh`:
```bash
sbatch bench_mpi.sh
```
This way, you will launch the code with different number of processes. 
This will write the runtimes of the strong scaling and the weak scaling in respectively `strong_scaling_mpi.csv` and `weak_scaling_mpi.csv`. You can then plot the results (strong scaling and weak scaling) with the python code `plot_mpi.py`