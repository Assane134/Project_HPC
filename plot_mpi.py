import numpy as np
import matplotlib.pyplot as plt

#Load and process strong scaling MPI data
strong = np.loadtxt('strong_scaling_mpi.csv', delimiter=',', skiprows=1)
procs_strong = strong[:, 0]
times_strong = strong[:, 1]
t1_strong = times_strong[0]
speedup_strong = t1_strong / times_strong
efficiency_strong = speedup_strong / procs_strong

# Plot strong scaling: Time vs Processes
plt.figure()
plt.plot(procs_strong, times_strong, 'o-', label='Temps d\'exécution (s)')
plt.xlabel('Nombre de processus')
plt.ylabel('Temps d\'exécution (secondes)')
plt.title('MPI Strong Scaling: Temps d\'exécution vs nombre de processus')
plt.grid(True)
plt.savefig('strong_scaling_mpi_time.png')
plt.close()


# Load and process weak scaling MPI data
weak = np.loadtxt('weak_scaling_mpi.csv', delimiter=',', skiprows=1)
procs_weak = weak[:, 0]
times_weak = weak[:, 2]
t1_weak = times_weak[0]
efficiency_weak = t1_weak / times_weak

# Plot weak scaling: Time vs Processes
plt.figure()
plt.plot(procs_weak, times_weak, 'o-', label='Temps d\'exécution (s)')
plt.xlabel('Nombre de processus')
plt.ylabel('Temps d\'exécution (secondes)')
plt.title('MPI Weak Scaling: Temps d\'exécution vs nombre de processus')
plt.grid(True)
plt.savefig('weak_scaling_mpi_time.png')
plt.close()
