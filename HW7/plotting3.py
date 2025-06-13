import matplotlib.pyplot as plt
import pandas as pd

# Read CSV files
openblas_data = pd.read_csv("openblas_performance.csv")
cublas_data = pd.read_csv("cublas_performance.csv")

# Merge data on 'n' to ensure alignment
data = pd.merge(openblas_data, cublas_data, on="n", how="outer")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(data["n"], data["OpenBLAS_GFLOPs"], label="OpenBLAS", marker="o", color="#1f77b4")
plt.plot(data["n"], data["cuBLAS_GFLOPs"], label="cuBLAS", marker="s", color="#ff7f0e")
plt.xlabel("Matrix Dimension (n)")
plt.ylabel("Performance (GFLOPs)")
plt.title("OpenBLAS vs cuBLAS Performance for Matrix Multiplication")
plt.xscale("log", base=2)
plt.grid(True, which="both", ls="--")
plt.legend()
plt.savefig("performance_plot.png", dpi=300, bbox_inches="tight")
plt.show()