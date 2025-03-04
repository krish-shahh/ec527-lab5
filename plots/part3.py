import matplotlib.pyplot as plt
import numpy as np

# Data from execution results
row_lengths = np.array([10, 12, 16, 22, 30, 40, 52, 66, 82, 100])
time_1_thread = np.array([26000, 26000, 48000, 92000, 168000, 296000, 502000, 814000, 1256000, 1870000])
time_2_threads = np.array([188000, 182000, 156000, 196000, 248000, 274000, 400000, 546000, 742000, 1118000])
time_4_threads = np.array([450000, 226000, 254000, 242000, 308000, 290000, 372000, 494000, 688000, 788000])

# Plot execution times
plt.figure(figsize=(8, 5))
plt.plot(row_lengths, time_1_thread, 'o-', label="1 Thread")
plt.plot(row_lengths, time_2_threads, 's-', label="2 Threads")
plt.plot(row_lengths, time_4_threads, 'x-', label="4 Threads")

plt.xlabel("Row Length (Array Size)")
plt.ylabel("Execution Time (Cycles)")
plt.title("Execution Time vs. Array Size for Different Thread Counts")
plt.legend()
plt.grid(True)
plt.show()
