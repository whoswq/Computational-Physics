import numpy as np
import matplotlib.pyplot as plt
from logistic_1 import *

steps = 5000
sp_num = 256
r_array = np.linspace(3.1, 4, 200)
x1_array = 0.5 * (1 + 1 / r_array + np.sqrt(1 - 2 / r_array - 3 /
                                            (r_array * r_array)))
x2_array = 0.5 * (1 + 1 / r_array - np.sqrt(1 - 2 / r_array - 3 /
                                            (r_array * r_array)))
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
plt.plot(r_array, x1_array, label="$x_1$", color="green")
plt.plot(r_array, x2_array, label="$x_2$", color="red")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="lightblue", s=0.05)
plt.legend()
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("4_1.pdf")
plt.show()