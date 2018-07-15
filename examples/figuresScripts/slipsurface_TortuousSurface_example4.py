from numpy import array
from pybimstab.slope import NaturalSlope
from pybimstab.bim import BlocksInMatrix
from pybimstab.slipsurface import TortuousSurface
terrainCoords = array(
    [[-2.49, 0.1, 1.7, 3.89, 5.9, 8.12, 9.87, 13.29, 20.29,
      21.43, 22.28, 23.48, 24.65, 25.17],
     [18.16, 17.88, 17.28, 15.73, 14.31, 13.58, 13, 3.61, 3.61,
      3.32, 2.71, 2.23, 1.21, 0.25]])
slope = NaturalSlope(terrainCoords)
bim = BlocksInMatrix(slopeCoords=slope.coords, blockProp=0.3,
                     tileSize=0.35, seed=123)
# Without a preferred path and smoothing the surface
surface = TortuousSurface(
    bim, dist1=5, dist2=15.78, heuristic='euclidean',
    reverseLeft=False, reverseUp=False, smoothFactor=2,
    preferredPath=None)
fig = surface.plot()
