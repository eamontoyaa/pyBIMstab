from numpy import array
from pybimstab.astar import Astar
grid = array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
              [1, 1, 0, 0, 1, 0, 1, 0, 0, 0],
              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
              [0, 0, 1, 1, 0, 1, 1, 1, 0, 0],
              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
              [0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
              [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
              [0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])
for heuristic in ['manhattan', 'euclidean']:
    astar = Astar(grid, startNode=(0, 0), goalNode=(9, 9),
                  heuristic=heuristic, reverseLeft=True,
                  reverseUp=True, preferredPath=None)
    fig = astar.plot()
