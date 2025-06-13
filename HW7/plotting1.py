import matplotlib.pyplot as plt

def parse_benchmark_file(filename):
    data = {"daxpy": [], "dgemv": [], "dgemm": []}
    current = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                if "daxpy" in line:
                    current = "daxpy"
                elif "dgemv" in line:
                    current = "dgemv"
                elif "dgemm" in line:
                    current = "dgemm"
            elif line:
                n, mflops = map(float, line.split())
                data[current].append((n, mflops))

    return data

def plot_performance(data):
    plt.figure(figsize=(10, 6))

    for label, points in data.items():
        ns, mflops = zip(*points)
        plt.plot(ns, mflops, label=label.upper(), linewidth=2)

    plt.xlabel("Problem Size n", fontsize=12)
    plt.ylabel("Performance (MFLOPs)", fontsize=12)
    plt.title("OpenBLAS Performance: L1 (DAXPY), L2 (DGEMV), L3 (DGEMM)", fontsize=14)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("openblas_performance_plot.png")
    plt.show()

if __name__ == "__main__":
    benchmark_data = parse_benchmark_file("openblas_perf.txt")
    plot_performance(benchmark_data)
