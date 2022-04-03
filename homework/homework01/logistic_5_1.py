import numpy as np
import matplotlib.pyplot as plt
from logistic_1 import *

steps = 10000
sp_num = 256
r_array = np.linspace(3.1, 3.58, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array,
                converged_points[:, i],
                color="blue",
                s=0.02,
                marker=".")
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("5_1.pdf")
plt.show()
