import numpy as np
import matplotlib.pyplot as plt
# plot trajectory adiabatic invariant, calculate berry phase
v = 0.001
traj = np.loadtxt("4_berryphase_v=%.3f.txt" % v)
t_array = traj[:, 0]
y_array = traj[:, 3]
p_array = traj[:, 4]
x_array = traj[:, 2]
g_array = traj[:, 1]
omega_array = traj[:, 6]
J_array = traj[:, 5]


def cal_berry_phase():
    phi_b = np.zeros(len(t_array))
    dt = t_array[1] - t_array[0]
    phi = 0
    for i in range(len(t_array)):
        phi += omega_array[i] * dt
        # for problem 4
        # phi_b[i] = (1 * 2.7556 - phi) % (2 * np.pi)
        # for problem 5
        # phi_b[i] = (1.000 * 3.6833 - phi) % (2 * np.pi)
        phi_b[i] = (1.000 * 1.347 - phi) % (2 * np.pi)
    print(phi_b[-1])


plt.plot(y_array, p_array, label="v = %.3f" % v, linewidth=0.2)
plt.xlabel("y")
plt.ylabel("p")
plt.title("trajectory")
plt.legend()
plt.savefig("4_traj_v=%.3f.pdf" % v)
plt.show()
plt.close()

plt.plot(x_array, J_array, label="v = %.3f" % v, linewidth=0.2)
plt.xlabel("g")
plt.ylabel("J")
plt.title("adiabatic invariant")
plt.legend()
plt.savefig("4x_ad_invr_v=%.3f.pdf" % v)
plt.show()
plt.close()

plt.plot(g_array, J_array, label="v = %.3f" % v, linewidth=0.2)
plt.xlabel("g")
plt.ylabel("J")
plt.title("adiabatic invariant")
plt.legend()
plt.savefig("4g_ad_invr_v=%.3f.pdf" % v)
plt.show()
plt.close()

cal_berry_phase()