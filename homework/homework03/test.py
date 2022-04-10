import numpy as np

N = 2000
gama = 2
alpha = 0.02
A = np.zeros((N, N))
eigvalue = np.sort(np.fromfile("test.dat", dtype="float64"))
for i in range(N):
    for j in range(i - 1, i + 2):
        if i == j:
            A[i, i] = gama * np.cos(np.pi * j * alpha * 2.0) + 2
        else:
            A[i, (j + N) % N] = -1
e, _ = np.linalg.eigh(A)
print(e)
print(eigvalue)
print(sum(abs(np.sort(eigvalue) - np.sort(e))))