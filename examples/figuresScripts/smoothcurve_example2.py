import numpy as np
from pybimstab.smoothcurve import SmoothCurve
x = np.linspace(0, 2*np.pi, 50)
y = np.sin(x) + np.random.random(50) * 0.5
for k in [2, 7, 15]:
    curve = SmoothCurve(x, y, k)
    fig = curve.plot()
