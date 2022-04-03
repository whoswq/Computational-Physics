import numpy as np
import matplotlib.pyplot as plt

r1 = np.linspace(0, 1, 50)
x1 = np.zeros(50)
r2 = np.linspace(1, 3, 100)
x2 = 1 - 1 / r2
plt.plot(r1, x1)
plt.plot(r2, x2)
plt.xlabel("$r$")
plt.ylabel("$x^*$")
plt.savefig("2.pdf")
plt.close()