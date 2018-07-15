from numpy import array
from pybimstab.slope import NaturalSlope
from pybimstab.watertable import WaterTable
from pybimstab.bim import BlocksInMatrix
from pybimstab.slipsurface import CircularSurface, TortuousSurface
from pybimstab.slices import MaterialParameters, Slices
terrainCoords = array(
    [[-2.49, 0.1, 1.7, 3.89, 5.9, 8.12, 9.87, 13.29, 20.29,
      21.43, 22.28, 23.48, 24.65, 25.17],
     [18.16, 17.88, 17.28, 15.73, 14.31, 13.58, 13, 3.61, 3.61,
      3.32, 2.71, 2.23, 1.21, 0.25]])
slope = NaturalSlope(terrainCoords)
bim = BlocksInMatrix(slopeCoords=slope.coords, blockProp=0.3,
                     tileSize=0.35, seed=123)
watertabDepths = array([[0, 5, 10, 15],
                        [8, 7, 3, 0]])
watertable = WaterTable(slopeCoords=slope.coords,
                        watertabDepths=watertabDepths,
                        smoothFactor=3)
preferredPath = CircularSurface(
    slopeCoords=slope.coords, dist1=5, dist2=15.78, radius=20)
surface = TortuousSurface(
    bim, dist1=4, dist2=15.78, heuristic='euclidean',
    reverseLeft=False, reverseUp=False, smoothFactor=2,
    preferredPath=preferredPath.coords, prefPathFact=2)
material = MaterialParameters(
    cohesion=15, frictAngle=23, unitWeight=17,
    blocksUnitWeight=21, wtUnitWeight=9.8)
slices = Slices(
    material=material, slipSurfCoords=surface.coords,
    slopeCoords=slope.coords, numSlices=10,
    watertabCoords=watertable.coords, bim=bim)
fig = slices.plot()
