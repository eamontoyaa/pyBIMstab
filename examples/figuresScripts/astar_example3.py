from numpy import array
from pybimstab.bim import BlocksInMatrix
from pybimstab.astar import Astar
seed = 111  # for repeatibilty
boundary = array([[-5, 0, 5, 0, -5], [0, 10, 0, -10, 0]])
bim = BlocksInMatrix(slopeCoords=boundary, blockProp=0.2,
                     tileSize=0.4, seed=seed)
astar = Astar(bim.grid, startNode=(0, 12), goalNode=(49, 12),
              heuristic='manhattan', reverseLeft=True,
              reverseUp=True, preferredPath=None)
fig = astar.plot()
