import pandas as pd
import matplotlib.pyplot as plt

# Load data
data1 = pd.read_csv("p4-4proc-cus.csv")
data2 = pd.read_csv("p4-4proc-mpi.csv")
data3 = pd.read_csv("p4-32proc-cus.csv")
data4 = pd.read_csv("p4-32proc-mpi.csv")

# Plot bandwidth vs message size
plt.figure(figsize=(10, 6))
plt.plot(data1.iloc[:, 0], data1.iloc[:, 1], color='b', marker='o', markersize=3, label='4 Processes, my_broadcast')
plt.plot(data2.iloc[:, 0], data2.iloc[:, 1], color='g', marker='o', markersize=3, label='4 Processes, MPI_Bcast')
plt.plot(data3.iloc[:, 0], data3.iloc[:, 1], color='k', marker='o', markersize=3, label='32 Processes, my_broadcast')
plt.plot(data4.iloc[:, 0], data4.iloc[:, 1], color='r', marker='o', markersize=3, label='32 Processes, MPI_Bcast')

# plt.xscale('log', base = 2)
plt.xlabel("Message Size (Bytes)")
plt.ylabel("Bandwidth (Bytes/Second)")
plt.title("MPI Broadcast Bandwidth vs. Message Size")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Problem4.png")
plt.show()
