import numpy as np
import matplotlib.pyplot as plt
from logistic_1 import *

steps = 100
num = 20
n_array = np.array([i for i in range(steps)])
r_array = [3.4, 3.5, 3.56, 3.567]
x0_array = np.linspace(0.1, 0.9, num)
x_array_2 = []  # 分别存储不同初值之下的序列集合
x_array_4 = []
x_array_8 = []
x_array_16 = []

for x0 in x0_array:
    x_array_2.append(generate_array(x0, r_array[0], steps))
    x_array_4.append(generate_array(x0, r_array[1], steps))
    x_array_8.append(generate_array(x0, r_array[2], steps))
    x_array_16.append(generate_array(x0, r_array[3], steps))

x_array_2 = np.array(x_array_2)
x_array_4 = np.array(x_array_4)
x_array_8 = np.array(x_array_8)
x_array_16 = np.array(x_array_16)

# 失败的尝试
"""x_average_2 = [sum(x_array_2[:,i])/num for i in range(steps)]
x_average_4 = [sum(x_array_4[:,i])/num for i in range(steps)]
x_average_8 = [sum(x_array_8[:,i])/num for i in range(steps)]
x_average_16 = [sum(x_array_16[:,i])/num for i in range(steps)]

plt.plot(n_array, x_average_2, label="r = 3.4")
plt.plot(n_array, x_average_4, label="r = 3.5")
plt.plot(n_array, x_average_8, label="r = 3.56")
plt.plot(n_array, x_average_16, label="r = 3.567")
plt.legend()
plt.xlabel("$n$")
plt.ylabel("$\\left<x\\right>$")
plt.show()"""

# 正确的作图
for i in range(len(x0_array)):
    plt.scatter(n_array, x_array_2[i], s=0.4)
plt.xlabel("n")
plt.ylabel("x")
plt.savefig("5_2.pdf")
plt.title("$r = 3.40$")
plt.show()
plt.close()

for i in range(len(x0_array)):
    plt.scatter(n_array, x_array_4[i], s=0.4)
plt.xlabel("n")
plt.ylabel("x")
plt.savefig("5_4.pdf")
plt.title("$r = 3.50$")
plt.show()
plt.close()

for i in range(len(x0_array)):
    plt.scatter(n_array, x_array_8[i], s=0.4)
plt.xlabel("n")
plt.ylabel("x")
plt.savefig("5_8.pdf")
plt.title("$r = 3.56$")
plt.show()
plt.close()

for i in range(len(x0_array)):
    plt.scatter(n_array, x_array_16[i], s=0.4)
plt.xlabel("n")
plt.ylabel("x")
plt.savefig("5_16.pdf")
plt.title("$r = 3.567$")
plt.show()
plt.close()
