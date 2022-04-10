import numpy as np
import matplotlib.pyplot as plt

eigvalue = np.sort(np.fromfile("5_eigenvalue_1200_alpha_10_11.dat", dtype="float64"))
n = len(eigvalue)
c = 2.0 - eigvalue

plt.scatter(np.arange(n),
            c,
            label="$\\alpha = 10/11$",
            marker=".",
            s=1,
            color="red")
plt.xlabel("index")
plt.ylabel("$2 - \\lambda$")

# plt.savefig("2_eigenvalue.pdf")
# n_array = np.arange((int(n / 2)))
# eig1 = (np.sqrt(4 + 4 * np.cos(2 * np.pi * n_array / n)**2))
# plt.scatter(n_array,
#             -(np.sort(-eig1)),
#             label="n=720 theoretical",
#             marker=".",
#             s=1,
#             color="blue")
# plt.scatter(n_array + int(n / 2),
#             -np.sort(eig1),
#             marker=".",
#             s=1,
#             color="blue")
# plt.savefig("3_compare.pdf")
plt.legend()
plt.savefig("5_alpha_10_11.pdf")
plt.show()
