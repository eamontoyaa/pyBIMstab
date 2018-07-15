import numpy as np
from pybimstab.polygon import Polygon
xC = [np.cos(theta) for theta in np.linspace(0, 2*np.pi, 100)]
yC = [np.sin(theta) for theta in np.linspace(0, 2*np.pi, 100)]
coords = np.array([xC, yC])
np.random.seed(123)
x = np.random.uniform(-1, 1, 100)
y = np.random.uniform(-1, 1, 100)
polygon = Polygon(coordinates=coords)
polygon.isinside(x=x, y=y, meshgrid=False, want2plot=True)
