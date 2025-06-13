import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("strassen.csv")

plt.plot(data.iloc[:,0], data.iloc[:,1], '-o', markersize = 3)
plt.title("Performance for Strassen Square Matrix Multiply (finding the next power of two before splitting and ending splits at n = 32)")
plt.xlabel("n dimensions (only even dimensions)")
plt.ylabel("FLOPs")
plt.grid()
plt.savefig("strassen_performance.png")
plt.show()