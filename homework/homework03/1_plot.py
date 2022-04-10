import numpy as np
import matplotlib.pyplot as plt

eigvalue = np.sort(np.fromfile("1_eigenvalue_1024.dat", dtype="float64"))
n = len(eigvalue)
n_array = np.arange(n)
eig_acu = np.sort(2 * (1 - np.cos(2 * np.pi * n_array / n)))

print(sum(abs(np.sort(eig_acu) - np.sort(eigvalue))))
plt.plot(eig_acu, label="analytical results")
plt.plot(eigvalue, label="numerical results")
plt.legend()
plt.xlabel("index of eigenvalue")
plt.ylabel("eigenvalue")
plt.savefig("1_compare_n=%d.pdf" % n)
plt.show()
