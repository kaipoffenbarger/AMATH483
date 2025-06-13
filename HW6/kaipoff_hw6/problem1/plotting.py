import pandas as pd
import matplotlib.pyplot as plt

partc = pd.read_csv("p1partc.csv")
partd = pd.read_csv("p1partd.csv")

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(8, 12))  # 2 rows, 1 column

# Plot 1 (part c)
axs[0].plot(partc.iloc[:, 0], partc.iloc[:, 1], color='b', marker='o', markersize=3)
axs[0].set_title("Time vs Thread Count in Computing an Integral")
axs[0].set_xlabel("Thread Count")
axs[0].set_ylabel("Time (seconds)")
axs[0].grid(True)

# Plot 2 (part d)
axs[1].plot(partd.iloc[:, 0], partd.iloc[:, 1], color='b', marker='o', markersize=3)
axs[1].set_title("Numerical Error vs. Number of Partition Points in Computing an Integral (1 thread)")
axs[1].set_xlabel("Number of Partition Points (log scale)")
axs[1].set_ylabel("Numerical Error (log scale)")
axs[1].set_xscale('log', base=10)
axs[1].set_yscale('log', base=10)
axs[1].grid(True)

# Adjust layout to increase vertical spacing between plots
plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Increased from default to add vertical space
plt.savefig('Problem1.png')
