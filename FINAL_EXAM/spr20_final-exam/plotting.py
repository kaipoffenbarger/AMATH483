import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv("grayscale_times.csv")

# Plot sequential and threaded times
plt.figure(figsize=(8, 5))
plt.plot(df["num_threads"], df["sequential_time_ms"], marker='o', label="Sequential Time")
plt.plot(df["num_threads"], df["threaded_time_ms"], marker='o', label="Threaded Time")

plt.xlabel("Number of Threads")
plt.ylabel("Execution Time (ms)")
plt.title("Grayscale Execution Time vs Number of Threads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("grayscale_execution_time.png")  # Save the plot as an image
plt.show()
