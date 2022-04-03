import numpy as np
import matplotlib.pyplot as plt

g_array_1 = np.linspace(0, 2, 2000)
x_array_1 = (1 - (0.5 * g_array_1)**(2 / 3))**(3 / 2)
g_array_2 = np.linspace(0, 1, 1000)
x_array_2 = (1 - (g_array_2)**(2 / 3))**(3 / 2)

plt.plot(g_array_1, x_array_1, label="$\\frac{m}{kl} = 0.5$")
plt.plot(g_array_2, x_array_2, label="$\\frac{m}{kl} = 1.0$")
plt.legend()
plt.xlabel("$g$")
plt.ylabel("$\\frac{x}{l}$")
plt.savefig("1_boundary.pdf")
plt.show()