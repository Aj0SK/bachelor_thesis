import math
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

n = 5
x1 = np.arange(0, 10, 0.001)
y = 2*np.arctan(n*x1)/math.pi
x = np.arange(0, 1, 0.0001)
diff = (2*np.arctan(n*(-10))/math.pi)
y1 = 2*np.arctan(n*(x1-10))/math.pi - diff

plt.plot(x, y)
plt.plot(x, y1)
plt.plot(x, np.linspace(0, -diff, 10000))

plt.xlabel("FPR")
plt.ylabel("TPR")

plt.axes().set_aspect('equal', adjustable='box')

plt.show()
