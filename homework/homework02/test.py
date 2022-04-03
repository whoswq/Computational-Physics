import numpy as np
import matplotlib.pyplot as plt

traj = np.loadtxt("test.txt")
x_array = traj[:, 1]
p_array = traj[:, 2]
plt.plot(x_array, p_array)
plt.show()