import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-1, 1, 0.05)
N = 1

while N <= 7:
	y = np.sin((N + 1 / 2) * np.pi * x) / (2 * np.sin(np.pi * x / 2))
	plt.plot(x, y)
	N = N + 1

plt.axvline(x=0, c="grey")
plt.axhline(y=0, c="grey")

plt.plot()
plt.xlabel("x")
plt.ylabel("delta")

plt.grid()
plt.show()