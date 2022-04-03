import numpy as np
import matplotlib.pyplot as plt
from logistic_1 import *


# 绘制周期一到周期二的分岔点
"""steps = 5000
sp_num = 128
r_array = np.linspace(2.997, 2.9992, 200)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_1.pdf")
plt.show()"""


# 绘制周期2到周期4的分岔点
"""steps = 5000
sp_num = 64
r_array = np.linspace(3.448, 3.4505, 200)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.1)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_2.pdf")
plt.show()"""


# 绘制周期4到周期8的分岔点
"""steps = 10000
sp_num = 128
r_array = np.linspace(3.54, 3.55, 200)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.05)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_3.pdf")
plt.show()"""


# 周期8到周期16的分岔点
"""steps = 10000
sp_num = 128
r_array = np.linspace(3.562, 3.568, 300)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.05)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_4.pdf")
plt.show()"""


# 周期16到周期32的分岔点
"""steps = 10000
sp_num = 256
r_array = np.linspace(3.5686, 3.5689, 300)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.01)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_5.pdf")
plt.show()"""


# 分岔现象的总体描述
steps = 10000
sp_num = 512
r_array = np.linspace(3.562, 3.629, 400)
converged_points = []

for r in r_array:
    converged_points.append(generate_array(0.5, r, steps)[-sp_num:])
converged_points = np.array(converged_points)
print("finished")
for i in range(sp_num):
    plt.scatter(r_array, converged_points[:, i], color="blue", s=0.002)
plt.xlabel("r")
plt.ylabel("x")
plt.savefig("6_another.pdf")
plt.show()