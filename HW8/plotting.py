import pandas as pd
import matplotlib.pyplot as plt

data_fftw = pd.read_csv("fftw_perf.csv")
data_cufft = pd.read_csv("cufft_perf.csv")

plt.figure(figsize=(10,6))
plt.plot(data_fftw["n"], data_fftw["flops"], label="FFTW", marker="o")
plt.plot(data_cufft["n"], data_cufft["flops"], label="CUFFT", marker="s")

# Set x-ticks to match the 'n' values so they align with data points
xticks = sorted(set(data_fftw["n"]).union(set(data_cufft["n"])))
plt.xticks(xticks)

plt.xlabel("Grid Dimension Size $n$ (real size is $n^3$, 3D plane is a cube)")
plt.ylabel("FLOPs")
plt.title("3D Gradient Performance on HYAK")
plt.legend()
plt.grid(True)
plt.savefig("gradient_perf_comparison.png")
plt.show()
