import numpy as np 
import matplotlib.pyplot as plt
from logistic_9_1 import *

steps = 5000
sp_num = 256
r_array = np.linspace(0.1, 1.75, 200)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_2.pdf")
plt.show()