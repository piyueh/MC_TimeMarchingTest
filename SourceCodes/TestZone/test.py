import numpy as np
import matplotlib.pyplot as plt

tmp = np.fromfile("test.txt", count=16002, dtype=np.float64, sep=" ")
T = tmp[0:16002:2]
G = tmp[1:16002:2]

plt.figure()
plt.plot(T, G, '.-')
plt.show()
