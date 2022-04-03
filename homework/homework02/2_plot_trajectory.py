import numpy as np
import matplotlib.pyplot as plt

v = 0.00390625
traj = np.loadtxt("3_trajectory_v=%.8f.txt" % v)
y_array = traj[:, 2]
p_array = traj[:, 3]
x_array = traj[:, 1]
J_array = traj[:, 4]
plt.plot(y_array, p_array, label="v = %.8f" % v)
plt.xlabel("y")
plt.ylabel("p")
plt.title("trajectory")
plt.legend()
plt.savefig("2_traj_v=%.8f.png" % v, dpi=800)
plt.close()

plt.plot(x_array, J_array, label="v = %.8f" % v)
plt.xlabel("x")
plt.ylabel("J")
plt.title("adiabatic invariant")
plt.legend()
plt.savefig("2_ad_invr_v=%.8f.png" % v, dpi=800)
plt.show()
