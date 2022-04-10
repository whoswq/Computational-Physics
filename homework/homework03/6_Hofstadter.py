import numpy as np
import matplotlib.pyplot as plt

M = 500
N = 720
eigvalue = np.fromfile("6_Hofstadter.dat", dtype="float64").reshape((M, N))
c = 2.0 - eigvalue
alpha = np.linspace(0, 1, M+1)

for i in range(M):
    plt.scatter(alpha[:M],
                c[:, i],
                marker=".",
                s=0.1,
                color="red")
plt.xlabel("alpha")
plt.ylabel("$2 - \\lambda$")
plt.savefig("6_Hofstadter_720_500.pdf")
plt.show()