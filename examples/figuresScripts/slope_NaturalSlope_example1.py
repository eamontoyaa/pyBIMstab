from numpy import array
from pybimstab.slope import NaturalSlope
terrainCoords = array(
    [[0, 10, 18, 28], [16.571, 16.571,  4.571,  4.571]])
slope = NaturalSlope(terrainCoords)
fig = slope.plot()
