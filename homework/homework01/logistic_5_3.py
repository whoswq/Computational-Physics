from hashlib import new
import numpy as np
import matplotlib.pyplot as plt
from logistic_1 import *

steps = 100
n_array = np.array([i for i in range(steps)])
r_array = [3.4, 3.5, 3.56, 3.567]
x0 = 0.1
x_array_2 = generate_array(x0, r_array[0], steps)
x_array_4 = generate_array(x0, r_array[1], steps)
x_array_8 = generate_array(x0, r_array[2], steps)
x_array_16 = generate_array(x0, r_array[3], steps)


def cesaro_sum(array):
    n = len(array)
    new_array = np.zeros(n, dtype="float64")
    for i in range(n):
        new_array[i] = sum(array[0:i + 1]) / (i + 1)
    return new_array

"""plt.plot(n_array, x_array_2, label="origin sequence")
plt.plot(n_array, cesaro_sum(x_array_2), label="cesaro sum")
plt.legend()
plt.xlabel("steps")
plt.ylabel("x")
plt.savefig("6_cesaro_p=2.pdf")
plt.show()"""

# 不同周期序列的切萨罗求和比较
plt.plot(n_array, cesaro_sum(x_array_2), label="r = 3.4")
plt.plot(n_array, cesaro_sum(x_array_4), label="r = 3.5")
plt.plot(n_array, cesaro_sum(x_array_8), label="r = 3.56")
plt.plot(n_array, cesaro_sum(x_array_16), label="r = 3.567")
plt.legend()
plt.xlabel("steps")
plt.ylabel("x")
plt.savefig("6_cesaro_compare_another.pdf")
plt.show()