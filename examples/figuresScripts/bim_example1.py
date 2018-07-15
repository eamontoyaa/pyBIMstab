from numpy import array
from pybimstab.bim import BlocksInMatrix
slopeCoords = array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0]])
bim = BlocksInMatrix(slopeCoords=slopeCoords, blockProp=0.5,
                     tileSize=0.1, seed=123)
fig = bim.plot()
