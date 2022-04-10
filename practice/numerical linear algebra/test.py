import numpy as np
N = 10
A = np.zeros((N, N))
for i in range(N):
    for j in range(i - 1, i + 2):
        if (i == j):
            A[i, i] = 2
        else:
            A[i, (j + N) % N] = -1
e, _ = np.linalg.eigh(A)
print(e)
