import numpy as np
import matplotlib.pyplot as plt

N = 40  # total steps


def logistic_map(x, r):
    return r * (2 - np.cosh(x))


def generate_array(x0, r, steps):
    l_array = np.zeros(steps, dtype="float64")
    for i in range(steps):
        l_array[i] = x0
        x0 = logistic_map(x0, r)
    return l_array

if __name__ == "__main__":
    x_array = np.array([i for i in range(N)])
    # select 4 different initial condition for r = 0.5, r = 1.5
    y_array_0_0 = generate_array(0.2, 0.5, N)
    y_array_0_1 = generate_array(0.4, 0.5, N)
    y_array_0_2 = generate_array(0.6, 0.5, N)
    y_array_0_3 = generate_array(0.8, 0.5, N)
    y_array_1_0 = generate_array(0.2, 1.5, N)
    y_array_1_1 = generate_array(0.4, 1.5, N)
    y_array_1_2 = generate_array(0.6, 1.5, N)
    y_array_1_3 = generate_array(0.8, 1.5, N)
    plt.plot(x_array, y_array_0_0, label="$x_0 = 0.2$")
    plt.plot(x_array, y_array_0_1, label="$x_0 = 0.4$")
    plt.plot(x_array, y_array_0_2, label="$x_0 = 0.6$")
    plt.plot(x_array, y_array_0_3, label="$x_0 = 0.8$")
    plt.legend()
    plt.xlabel("steps")
    plt.ylabel("x")
    plt.title("$r = 0.5$")
    plt.savefig("9_1_0.pdf")
    plt.close()
    plt.plot(x_array, y_array_1_0, label="$x_0 = 0.2$")
    plt.plot(x_array, y_array_1_1, label="$x_0 = 0.4$")
    plt.plot(x_array, y_array_1_2, label="$x_0 = 0.6$")
    plt.plot(x_array, y_array_1_3, label="$x_0 = 0.8$")
    plt.legend()
    plt.xlabel("steps")
    plt.ylabel("x")
    plt.title("$r = 1.5$")
    plt.savefig("9_1_1.pdf")
    plt.show()
