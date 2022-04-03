import numpy as np 
import matplotlib.pyplot as plt
from logistic_9_1 import *


# 周期1到周期2
"""steps = 5000
sp_num = 64
r_array = np.linspace(1.152, 1.158, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_1.pdf")
plt.show()"""


# 周期2到周期4
"""steps = 5000
sp_num = 64
r_array = np.linspace(1.45, 1.455, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_2.pdf")
plt.show()"""

# 周期4到周期8
"""steps = 5000
sp_num = 128
r_array = np.linspace(1.50, 1.52, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_3.pdf")
plt.show()"""


# 周期8到周期16
"""steps = 10000
sp_num = 256
r_array = np.linspace(1.52, 1.526, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_4.pdf")
plt.show()"""


# 周期16到周期32
"""steps = 10000
sp_num = 256
r_array = np.linspace(1.5259, 1.5261, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_5.pdf")
plt.show()"""


steps = 10000
sp_num = 512
r_array = np.linspace(1.526, 1.54, 300)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.005)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("9_3_another.pdf")
plt.show()