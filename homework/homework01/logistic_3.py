from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import logistic_1

N = 20

r = 3.1

x_array = np.array([i for i in range(N)])
y_array_0 = logistic_1.generate_array(0.2, r, N)
y_array_1 = logistic_1.generate_array(0.4, r, N)
y_array_2 = logistic_1.generate_array(0.6, r, N)
y_array_3 = logistic_1.generate_array(0.8, r, N)
y_array_4 = logistic_1.generate_array(1.1, r, N)
plt.plot(x_array, y_array_0, label="$x_0 = 0.2$")
plt.plot(x_array, y_array_1, label="$x_0 = 0.4$")
plt.plot(x_array, y_array_2, label="$x_0 = 0.6$")
plt.plot(x_array, y_array_3, label="$x_0 = 0.8$")
plt.xlabel("steps")
plt.ylabel("x")
plt.title("$r = 3.1$")
plt.legend()
plt.savefig("3.pdf")
plt.show()