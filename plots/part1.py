import numpy as np
import matplotlib.pyplot as plt

# Data extracted from user's provided output
omega_values = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]

iterations_32x32 = [346, 274, 224, 187, 161, 142, 125, 113, 102, 94]
iterations_64x64 = [1212, 949, 777, 645, 551, 478, 428, 386, 350, 323]
iterations_128x128 = [4471, 3577, 2930, 2462, 2134, 1862, 1661, 1484, 1342, 1227]
iterations_256x256 = [17265, 14008, 11558, 9734, 8373, 7245, 6413, 5733, 5191, 4737]

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(omega_values, iterations_32x32, marker='o', linestyle='-', label="32x32 Grid")
plt.plot(omega_values, iterations_64x64, marker='s', linestyle='-', label="64x64 Grid")
plt.plot(omega_values, iterations_128x128, marker='^', linestyle='-', label="128x128 Grid")
plt.plot(omega_values, iterations_256x256, marker='d', linestyle='-', label="256x256 Grid")

# Graph aesthetics
plt.xlabel("Relaxation Parameter (ω)")
plt.ylabel("Iterations to Convergence")
plt.title("SOR Iterations vs. Relaxation Parameter for Different Grid Sizes")
plt.legend()
plt.grid(True)
plt.show()
