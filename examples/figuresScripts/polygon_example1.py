from numpy import array
from pybimstab.polygon import Polygon
coords = array([[0, 1, 1, 0], [0, 0, 1.5, 1.5]])
x, y = 0.5, 2
polygon = Polygon(coordinates=coords)
polygon.isinside(x=x, y=y, meshgrid=False, want2plot=True)
