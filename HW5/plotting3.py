import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_row = pd.read_csv("row_swap.csv")
data_col = pd.read_csv("col_swap.csv")

plt.plot(data_row.iloc[:,0], data_row.iloc[:,1], '-o', markersize = 3, color = 'b', label = "Row Swap")
plt.plot(data_col.iloc[:,0], data_col.iloc[:,1], '-o', markersize = 3, color = 'r', label = "Column Swap")
plt.title("Row Swap and Column Swap time on Column-Major Square Matrix")
plt.xlabel("n dimensions (powers of 2 starting at 16, ending at 4096)")
plt.ylabel("Time (seconds)")
plt.yscale("log", base = 10)
plt.grid()
plt.legend()
plt.savefig("swap_performance.png")
plt.show()