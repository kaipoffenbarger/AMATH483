import matplotlib.pyplot as plt
import pandas as pd

# Load data
data1 = pd.read_csv("FLOPS_1.csv")
data2 = pd.read_csv("FLOPS_2.csv")
data3 = pd.read_csv("FLOPS_3.csv")

# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(8, 12))  # 3 rows, 1 column

# Plot 1
axs[0].plot(data1.iloc[:, 0], data1.iloc[:, 1], marker='o')
axs[0].set_title("Performance Plot - Level 1 BLAS")
axs[0].set_xlabel(data1.columns[0])
axs[0].set_ylabel("FLOPs")
axs[0].grid(True)

# Plot 2
axs[1].plot(data2.iloc[:, 0], data2.iloc[:, 1], marker='s', color='orange')
axs[1].set_title("Performance Plot - Level 2 BLAS")
axs[1].set_xlabel(data2.columns[0])
axs[1].set_ylabel("FLOPs")
axs[1].grid(True)

# Plot 3
axs[2].plot(data3.iloc[:, 0], data3.iloc[:, 1], marker='^', color='green')
axs[2].set_title("Performance Plot - Level 3 BLAS")
axs[2].set_xlabel(data3.columns[0])
axs[2].set_ylabel("FLOPs")
axs[2].grid(True)

plt.tight_layout()
plt.show()

