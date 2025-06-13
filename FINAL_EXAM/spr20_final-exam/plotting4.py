import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV data
df = pd.read_csv("output_times.csv")

# Compute average wall time for each process count
avg_times = df.groupby("processes")["wall_time"].mean().reset_index()

# Extract data
processes = avg_times["processes"].values
wall_times = avg_times["wall_time"].values

# Ideal scaling based on time at 1 process
ideal_time = wall_times[0]
ideal_scaling = ideal_time / processes

# Plot
plt.figure(figsize=(8, 5))
plt.plot(processes, wall_times, marker='o', label="Measured Wall Time")
plt.plot(processes, ideal_scaling, linestyle='--', color='gray', label="Ideal Scaling")

plt.xlabel("Number of Processes")
plt.ylabel("Wall Time (seconds)")
plt.title("Strong Scaling of Parallel Ax = b Solver (N=4096)")
plt.xticks(processes)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("solver_scaling_plot.png")  # Save the plot as an image
plt.show()
