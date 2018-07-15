from numpy import array
from pybimstab.astar import PreferredPath, Astar
grid = array([[0,  0,  0,  0,  1,  0,  0,  0,  0,  0],
              [0,  0,  0,  0,  1,  1,  0,  0,  0,  0],
              [0,  0,  0,  1,  1,  0,  1,  0,  1,  0],
              [0,  1,  1,  0,  0,  0,  0,  1,  0,  0],
              [0,  0,  0,  1,  0,  0,  0,  1,  1,  1],
              [0,  0,  1,  0,  1,  1,  0,  0,  0,  0],
              [0,  0,  1,  0,  1,  0,  1,  0,  0,  0],
              [1,  1,  0,  1,  0,  0,  0,  0,  1,  1],
              [1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
              [0,  0,  0,  0,  1,  0,  1,  0,  0,  1],
              [0,  1,  0,  0,  0,  1,  0,  0,  0,  0],
              [1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
              [0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
              [1,  0,  0,  0,  1,  0,  0,  0,  0,  0]
              ])
coordsIdx = array(
    [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 13,
      13, 13, 13, 13, 13, 13, 13],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4,
      5, 6, 7, 8, 9]])
preferredPath = PreferredPath(coordsIdx, factor=1)
for typePP in [None, preferredPath]:
    astar = Astar(grid, startNode=(0, 0), goalNode=(13, 9),
                  heuristic='manhattan', reverseLeft=True,
                  reverseUp=True, preferredPath=typePP)
    fig = astar.plot(plotPreferredPath=True)
