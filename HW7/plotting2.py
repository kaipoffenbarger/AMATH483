import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt('zgesv_errors.txt', comments='#')
n = data[:, 0]
log_residual = data[:, 1]
log_normalized_error = data[:, 2]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(n, log_residual, 'o-', label='log10(residual)', linewidth=2)
plt.plot(n, log_normalized_error, 's--', label='log10(normalized error)', linewidth=2)

plt.xlabel('Matrix size n')
plt.ylabel('log10(Error)')
plt.title('LAPACKE_zgesv Accuracy vs Matrix Size')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig('zgesv_error_plot.png', dpi=300)
plt.show()
