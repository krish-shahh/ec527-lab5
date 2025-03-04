import numpy as np
import matplotlib.pyplot as plt

# Data from the output
grid_sizes = np.array([34, 58, 98, 154, 226])
sor_times = np.array([8.14e5, 9.212e6, 5.381e7, 3.615e8, 7.046e8])  # Standard SOR time (ns)
sor_iters = np.array([62, 210, 385, 1033, 908])

redblack_times = np.array([4.62e5, 3.67e6, 3.116e7, 1.902e8, 2.532e8])  # Red/Black SOR time (ns)
redblack_iters = np.array([93, 249, 637, 672, 1015])

reversed_times = np.array([1.422e6, 1.542e7, 8.786e7, 5.942e8, 1.149e9])  # Reversed SOR time (ns)
reversed_iters = np.array([62, 210, 385, 1033, 908])

blocked_times = np.array([5.52e5, 5.614e6, 3.124e7, 2.13e8, 4.005e8])  # Blocked SOR time (ns)
blocked_iters = np.array([62, 210, 385, 1033, 908])

# Compute time per innermost loop iteration
sor_time_per_iter = sor_times / (sor_iters * grid_sizes**2)
redblack_time_per_iter = redblack_times / (redblack_iters * grid_sizes**2)
reversed_time_per_iter = reversed_times / (reversed_iters * grid_sizes**2)
blocked_time_per_iter = blocked_times / (blocked_iters * grid_sizes**2)

# Plot the results
plt.figure(figsize=(8, 5))
plt.plot(grid_sizes, sor_time_per_iter, marker='o', linestyle='-', label="Standard SOR")
plt.plot(grid_sizes, redblack_time_per_iter, marker='s', linestyle='--', label="Red/Black SOR")
plt.plot(grid_sizes, reversed_time_per_iter, marker='d', linestyle='-.', label="Reversed SOR")
plt.plot(grid_sizes, blocked_time_per_iter, marker='x', linestyle=':', label="Blocked SOR")

plt.xlabel("Grid Size")
plt.ylabel("Time per Iteration (ns)")
plt.title("Time per Innermost Loop Iteration for SOR Methods")
plt.legend()
plt.grid()
plt.show()
