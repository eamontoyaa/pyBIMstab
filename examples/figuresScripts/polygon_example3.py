from numpy import array
from pybimstab.polygon import Polygon
coords = array([[0, 1, 1, 0.], [0, 0, 1.5, 1.5]])
x = array([0.3, 0.5, 0.7, 1.2, 1.0])
y = array([0.6, 1., 1.4, 2.4, 0])
polygon = Polygon(coordinates=coords)
polygon.isinside(x=x, y=y, meshgrid=True, want2plot=True)
