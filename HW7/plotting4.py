import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
try:
    data = pd.read_csv("bandwidth.csv")
except FileNotFoundError:
    print("Error: bandwidth.csv not found. Ensure the file exists in the current directory.")
    exit(1)

# Verify expected columns
expected_columns = ["BufferSize_Bytes", "H2D_Bandwidth_GBs", "D2H_Bandwidth_GBs"]
if not all(col in data.columns for col in expected_columns):
    print(f"Error: bandwidth.csv must contain columns: {expected_columns}")
    exit(1)


# Plot
plt.figure(figsize=(10, 6))
plt.plot(data["BufferSize_Bytes"], data["H2D_Bandwidth_GBs"], label="Host-to-Device", marker="o", color="#1f77b4")
plt.plot(data["BufferSize_Bytes"], data["D2H_Bandwidth_GBs"], label="Device-to-Host", marker="s", color="#ff7f0e")
plt.xlabel("Buffer Size in Bytes")
plt.ylabel("Bandwidth (GB/s)")
plt.title("Host-to-Device vs Device-to-Host Bandwidth for Data Copies")
plt.xscale("log", base=2)
plt.grid(True, which="both", ls="--")
plt.legend()
plt.savefig("bandwidth_plot.png", dpi=300, bbox_inches="tight")
plt.show()