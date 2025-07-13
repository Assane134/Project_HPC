import numpy as np
import matplotlib.pyplot as plt

# Load benchmarking data
data = np.loadtxt('openmp_timings.csv', delimiter=',', skiprows=1)
threads = data[:, 0]
times   = data[:, 1]

# Compute metrics
t_seq      = times[threads == 1][0]
speedup    = t_seq / times
efficiency = speedup / threads

# Plot Speedup and save
plt.figure(figsize=(6, 4))
plt.plot(threads, speedup, 'o-', label='Speedup')
plt.xlabel('Nombre de threads')
plt.ylabel('Speedup')
plt.title('Speedup vs Threads OpenMP')
plt.grid(True)
plt.tight_layout()
plt.savefig('speedup_plot.png', dpi=300)
plt.close()

# Plot Efficiency and save
plt.figure(figsize=(6, 4))
plt.plot(threads, efficiency, 's--', label='Efficiency', color='C1')
plt.xlabel('Nombre de threads')
plt.ylabel('Efficiency')
plt.ylim(0.2, 1.05)
plt.title('Efficiency vs Threads OpenMP')
plt.grid(True)
plt.tight_layout()
plt.savefig('efficiency_plot.png', dpi=300)
plt.close()
