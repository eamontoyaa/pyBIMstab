# -*- coding: utf-8 -*-
'''
Module for defining the classes related to the :math:`\\mathrm{A}^\\ast`
algorithm.

:math:`\\mathrm{A}^\\ast` algorithm (pronounced as "A star") is a search
algorithm proposed in 1968 by
`Hart, Nilsson & Raphael (1968) <https://doi.org/10.1109/TSSC.1968.300136>`_.
It optimally finds the minimum cost path in a graph from a start node
:math:`s` to a goal node :math:`g`; in other words, the shortest path if the
cost is given by its length. In addition, the algorithm does not expand as many
nodes as other algorithms do; then, for a graph with a huge quantity of nodes,
the computational performance is higher with respect to other search algorithms

The :math:`\\mathrm{A}^\\ast` algorithm works constantly evaluating an
evaluation function :math:`f(n)` composed of two parts: one is the actual cost
of the optimum path traveled from :math:`s` to the current node :math:`n`
given by the expression :math:`g(n)`, and the second is the cost of the optimum
path from the current node :math:`n` to :math:`g` given by the expression
:math:`h(n)`, which is the heuristic component of the algorithm and it could be
for example either the
`euclidean <https://en.wikipedia.org/wiki/euclidean_distance>`_ or
`manhattan <https://en.wikipedia.org/wiki/Taxicab_geometry>`_ distance.
Thereby, the evaluation function to measure the path cost is
:math:`f(n) = g(n) + h(n)`.
'''


# %%
class PreferredPath:
    '''Creates an instance of an object that stores the indexes of a
    polyline in the space of the grid-graph that represents a preferential path
    to be followed when the :math:`\\mathrm{A}^\\ast` algorithm is applied. ::

        PreferredPath(coordsIdx, factor)

    The object has an attribute called factor that represents a coefficient
    :math:`k` that multiplies the distance :math:`d` between the current node
    :math:`n` and the polyline. Considering the above, the function for
    evaluating the total cost of a node is modified as
    :math:`f(n) = g(n) + h(n) + kd`

    Attributes:
        polyline (`numpy.ndarray`): (2, n) array with the indexes of a polyline
            in the space of a grid-graph where the :math:`\\mathrm{A}^\\ast`
            algorithm is applied. The first row corresponds to the rows and the
            second to the columns of the grid-graph.
        factor (`int` or `float`): Multiplier of the shortest distance
            between the current node and the polyline.

    Examples:
        >>> from numpy import array
        >>> from pybimstab.astar import PreferredPath
        >>> coordsIdx = array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                13, 13, 13, 13, 13, 13, 13, 13, 13, 13],
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                2, 3, 4, 5, 6, 7, 8, 9]])
        >>> preferredPath = PreferredPath(coordsIdx, factor=1)
        >>> preferredPath.__dict__
        {'coordsIdx': array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13,
                              13, 13, 13, 13, 13, 13, 13, 13, 13],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                              2, 3, 4, 5, 6, 7, 8, 9]]),
         'factor': 1}
    '''
    def __init__(self, coordsIdx, factor):
        '''
        PreferredPath(coordsIdx, factor)
        '''
        from numpy import array

        self.coordsIdx = array(coordsIdx)
        self.factor = factor


# %%
class Node:
    '''Creates an instance of an object that defines the structure of a
    node that belongs to a grid-graph where is wanted to find the optimum path
    through the :math:`\\mathrm{A}^\\ast` algorithm. ::

        Node(pos=None, father=None, gCost=None, hCost=None, val=0)

    Attributes:
        pos (`tuple`): Indexes of the the current node in the grid-graph.
        father (`tuple`): Indexes of the father node of the current node.
            Default value is ``None``.
        gCost (`int` or `float`): Length of the traveled path from the
            start node to the current node. Default value is ``None``.
        hCost (`int` or `float`): Heuristic length of the path from the
            current node to the goal node. Default value is ``None``.
        val (`int`): Value that store the node cell in the matrix that defines
            the grid-graph. Default value is ``0``.

    Note:
        The class ``Node`` requires `numpy <http://www.numpy.org/>`_ and
        `shapely <https://pypi.python.org/pypi/Shapely>`_

    Examples:
        >>> node = Node(pos=(6, 7), father=(5, 7), gCost=5, hCost=10, val=1)
        >>> node.__dict__
        {'father': (5, 7), 'gCost': 5, 'hCost': 10, 'pos': (6, 7), 'val': 1}
        '''
    def __init__(self, pos, father=None, gCost=None, hCost=None, val=0):
        '''
        Node(pos=None, father=None, gCost=None, hCost=None, val=0)
        '''
        self.pos = pos
        self.father = father
        self.gCost = gCost
        self.hCost = hCost
        self.val = val

    def getHcost(self, goalNode, heuristic='manhattan', preferredPath=None):
        '''Method to obtain the heuristic component, :math:`h(n)`, that
        estimates the cost (or length) of the shortest path from the current
        node :math:`n` to the goal node.

        It must be selected either ``manhattan`` or ``euclidean`` as the model
        to estimate the length of the optimum path.

            - **manhattan** is the sum of the cathetus of the right triangle\
                defined by the current node and the goal node.
            - **euclidean** is the length of the hypotenuse of the right\
                triangle defined by the current node and the goal node.

        It is possible to append a polyline to force the path to follow a
        preferential path.

        Args:
            goalNode (`Node` object): object with the structure of the goal
                node.
            heuristic (`str`): Name of the geometric model to determine the
                heuristic distance. It must be selected either ``manhattan``
                or ``euclidean``. The first one is the default value.
            preferredPath (`Node` object): Optional argument of the class
                ``Node`` to force the path. ``None`` is the default value.

        Returns:
            (`int` or `float`): value of the estimated heuristic distance of\
                the opmtimum path.

        Examples:
            >>> from pybimstab.astar import Node
            >>> goalNode = Node(pos=(9, 9))
            >>> node = Node(pos=(0, 0))
            >>> node.getHcost(goalNode, heuristic='manhattan',
            >>>               preferredPath=None)
            18
            >>> from pybimstab.astar import Node
            >>> goalNode = Node(pos=(9, 9))
            >>> node = Node(pos=(0, 0))
            >>> node.getHcost(goalNode, heuristic='euclidean',
            >>>               preferredPath=None)
            12.727922061357855
        '''
        from shapely.geometry import LineString, Point

        # Verifying if is wanted to append the extra modifier to the heuristic
        plusCost = 0
        if preferredPath is not None:
            pline = LineString(preferredPath.coordsIdx.T)
            point = Point(self.pos)
            plusCost = point.distance(pline) * preferredPath.factor
        # extract coordinates from structure
        cNode = self.pos
        gNode = goalNode.pos
        if heuristic == 'manhattan':
            cost = abs(cNode[0]-gNode[0]) + abs(cNode[1]-gNode[1])
        elif heuristic == 'euclidean':
            cost = ((cNode[0]-gNode[0])**2 + (cNode[1]-gNode[1])**2)**0.5
        else:
            print('Invalid heuristic, You must have selected manhattan or ' +
                  'euclidean. manhattan selected by default')
            cost = abs(cNode[0]-gNode[0]) + abs(cNode[1]-gNode[1])
        cost += plusCost
        setattr(self, 'hCost', cost)
        return cost

    def getGcost(self, fatherNode):
        '''Method to obtain the cumulated cost :math:`g(n)`, of the traveled
        path from the start node to the current node :math:`n`.

        Args:
            fatherNode (``Node`` object): object with the structure\
                of the current node's father.

        Returns:
            (`int` or `float`): traveled-path length from the the start node\
            to the current node.

        Examples:
            >>> from pybimstab.astar import Node
            >>> node = Node(pos=(9, 9))
            >>> fatherNode = Node(pos=(9, 8), gCost=15)
            >>> node.getGcost(fatherNode)
            16

            >>> from pybimstab.astar import Node
            >>> node = Node(pos=(9, 9))
            >>> fatherNode = Node(pos=(8, 8), gCost=15)
            >>> node.getGcost(fatherNode)
            16.4142
        '''
        from numpy import array

        if any(array(self.pos) == array(fatherNode.pos)):
            cost = fatherNode.gCost + 1  # horizontal or vertical movement.
        else:
            cost = fatherNode.gCost + 1.4142  # diagonal movement
        setattr(self, 'gCost', cost)
        return cost


# %%
class Astar:
    '''Creates an instance of an object that defines the optimum path into a
    grid-graph maze, from a start node to a goal node given. ::

        Astar(grid, startNode, goalNode, heuristic='manhattan',
              reverseLeft=True, reverseUp=True, preferredPath=None)

    Attributes:
        gridGraph (`MazeStructure` object): object with the structure\
            of maze where is wanted to find the optimum path.
        startNode (`tuple`, `list` or `numpy.ndarray`): indexes of \
            the ``gridGraph.matrix`` where the initial node is located. It has\
            to be a matrix cell, *ie*, ``gridGraph.matrix[startNode]==0``.
        goalNode (`tuple`, `list` or `numpy.ndarray`): indexes of the\
            ``gridGraph.matrix`` where the ending node is located. It has to\
            be a matrix cell, *ie*, ``gridGraph.matrix[goalNode]==0``.
        heuristic (`str`): Name of the geometric model to determine the\
                heuristic distance. It must be selected either ``manhattan``\
                or ``euclidean``. The first one is the default value.
        reverseLeft (`bool`): Logical variable to allow or not reverses\
            movements to the left. Default value is ``True``.
        reverseUp (`bool`): Logical variable to allow or not reverses\
            movements to upward. Default value is ``True``.
        *forcedPath (`PreferredPath` object): Optional aguments to force the\
            optimum path close to a specific polyline.

    Note:
        The class ``Astar`` requires `NumPy <http://www.numpy.org/>`_\
        and `Matplotlib <https://matplotlib.org/>`_.

    Examples:
        >>> from numpy import array
        >>> from pybimstab.astar import PreferredPath, Astar
        >>> grid = array([[ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0],
        >>>               [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0],
        >>>               [ 0,  0,  0,  1,  1,  0,  1,  0,  1,  0],
        >>>               [ 0,  1,  1,  0,  0,  0,  0,  1,  0,  0],
        >>>               [ 0,  0,  0,  1,  0,  0,  0,  1,  1,  1],
        >>>               [ 0,  0,  1,  0,  1,  1,  0,  0,  0,  0],
        >>>               [ 0,  0,  1,  0,  1,  0,  1,  0,  0,  0],
        >>>               [ 1,  1,  0,  1,  0,  0,  0,  0,  1,  1],
        >>>               [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        >>>               [ 0,  0,  0,  0,  1,  0,  1,  0,  0,  1],
        >>>               [ 0,  1,  0,  0,  0,  1,  0,  0,  0,  0],
        >>>               [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        >>>               [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
        >>>               [ 1,  0,  0,  0,  1,  0,  0,  0,  0,  0]
        >>>               ])
        >>> coordsIdx = array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        >>>                     13, 13, 13, 13, 13, 13, 13, 13, 13, 13],
        >>>                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        >>>                     2, 3, 4, 5, 6, 7, 8, 9]])
        >>> preferredPath = PreferredPath(coordsIdx, factor=1)
        >>> astar = Astar(grid, startNode=(0, 0), goalNode=(13, 9),
        >>>               heuristic='manhattan', reverseLeft=True,
        >>>               reverseUp=True, preferredPath=preferredPath)
        >>> astar.__dict__.keys()
        dict_keys(['grid', 'startNode', 'goalNode', 'heuristic', 'reverseLeft',
                   'reverseUp', 'preferredPath', 'mazeStr', 'optimumPath'])
        '''

    def __init__(self, grid, startNode, goalNode, heuristic='manhattan',
                 reverseLeft=True, reverseUp=True, preferredPath=None):
        '''
        Astar(grid, startNode, goalNode, heuristic='manhattan',
              reverseLeft=True, reverseUp=True, preferredPath=None)
        '''
        self.grid = grid
        self.startNode = startNode
        self.goalNode = goalNode
        self.heuristic = heuristic
        self.reverseLeft = reverseLeft
        self.reverseUp = reverseUp
        self.preferredPath = preferredPath
        # for defining the maze structure
        self.defineMazeStr()
        # Get the optimum path via A star algortihm
        self.getPath()

    def defineMazeStr(self):
        '''Defines the clean structure of the grid-graph using objects
        instaced from the class ``Node`` for each node (or cell of the grid).

        Those cells with value equal to one or minus one will have a g-value
        equal to infinite because over those cells, a path is impossible to be
        traced.

        Returns:
            (`numpy.ndarray`): array with the same shape of the input grid\
                where each cell is an object from the class ``Node``\
                with the structure of its respective node.

        Examples:
            >>> # after executing the example of the class
            >>> astar.mazeStr[0, 0].__dict__
            {'father': (0, 0), 'gCost': 0, 'hCost': 22.0, 'pos': (0, 0),
             'val': 0}
        '''
        import numpy as np

        numberRows, numberCols = self.grid.shape
        mazeStr = np.empty((numberRows, numberCols), dtype=object)
        for m in range(numberRows):
            for n in range(numberCols):
                mazeStr[m, n] = Node(pos=(m, n), val=self.grid[m, n])
                # blocks or outside cells have infinite value
                if abs(self.grid[m, n]):
                    mazeStr[m, n].gCost = np.inf
        setattr(self, 'mazeStr', mazeStr)

        # setting the structure to the start node
        mazeStr[self.startNode].father = self.startNode
        mazeStr[self.startNode].getHcost(goalNode=mazeStr[self.goalNode],
                                         heuristic=self.heuristic,
                                         preferredPath=self.preferredPath)
        mazeStr[self.startNode].gCost = 0
        return mazeStr

    def __checkerrors(self):
        '''Function to check if there is any problem with either the start
        or goal nodes. It includes either the existence of a block or an
        unallowed cell in those nodes.'''
        if self.mazeStr[self.startNode].val == -1:
            raise ValueError('Start node: not in the allowed cells')
        elif self.mazeStr[self.startNode].val == 1:
            raise ValueError('Start node: is a block-cell')
        elif self.mazeStr[self.goalNode].val == -1:
            raise ValueError('Goal node: not in the allowed cells')
        elif self.mazeStr[self.goalNode].val == 1:
            raise ValueError('Goal node: is a block-cell')
        return

    def getNeighbours(self, node):
        '''Method for obtaining the possible neighbours of an specific node
        that belongs to the grid given.

        Each neighbour is given as a tuple with the indexes of the grid.

        Args:
            node (``Node`` object): object with the structure of the node which
                is wanted to know its possible neighbours.

        Returns:
            (`list`): Tuples with the indexes of possible neighbours of the\
                node in question.

        Examples:
            >>> # after executing the example of the class
            >>> from pybimstab.astar import Node
            >>> astar.getNeighbours(Node(pos=(1, 1)))
            [(0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0), (0, 0)]
            >>> astar.getNeighbours(Node(pos=(0, 0)))
            [(0, 1), (1, 1), (1, 0)]
            >>> astar.getNeighbours(Node(pos=(0, 1)))
            [(0, 2), (1, 2), (1, 1), (1, 0), (0, 0)]
            >>> astar.getNeighbours(Node(pos=(0, 2)))
            [(1, 2), (1, 1), (0, 1)]
            >>> astar.getNeighbours(Node(pos=(1, 2)))
            [(0, 2), (2, 2), (2, 1), (1, 1), (0, 1)]
            >>> astar.getNeighbours(Node(pos=(2, 2)))
            [(1, 2), (2, 1), (1, 1)]
            >>> astar.getNeighbours(Node(pos=(2, 1)))
            [(1, 1), (1, 2), (2, 2), (2, 0), (1, 0)]
            >>> astar.getNeighbours(Node(pos=(2, 0)))
            [(1, 0), (1, 1), (2, 1)]
            >>> astar.getNeighbours(Node(pos=(1, 0)))
            [(0, 0), (0, 1), (1, 1), (2, 1), (2, 0)]
        '''
        i, j = node.pos
        numberRows, numberCols = self.grid.shape
        allNeighbours = [(i-1, j), (i-1, j+1), (i, j+1), (i+1, j+1), (i+1, j),
                         (i+1, j-1), (i, j-1), (i-1, j-1)]
        neighbours = list()
        for neighbour in allNeighbours:
            if neighbour[0] < 0 or neighbour[1] < 0 or \
                    neighbour[0] == numberRows or neighbour[1] == numberCols:
                continue
            if not self.reverseLeft and neighbour[1] < j:
                continue
            if not self.reverseUp and neighbour[0] > i:
                continue
            neighbours.append(neighbour)
        return neighbours

    def getWayBack(self, node):
        '''Method for obtaining the whole way back of an specific node which
        has been opened by the *A star* algorithm.

        Args:
            node (``Node`` object): object with the structure of the node which
                is wanted to know its way back.

        Returns:
            (`numpy.ndarray`): :math:`\\left(2 \\times n \\right)` array\
                where :math:`n` is the number of nodes where the path has\
                crossed; the first row of the array contains the abscisses and\
                the second one contains the ordinates of the nodes into the\
                grid-graph

        Examples:
            >>> import numpy as np
            >>> from pybimstab.astar import PreferredPath, Node, Astar
            >>> grid = np.zeros((3,3))
            >>> grid[1:3, 0:2] = 1
            >>> astar = Astar(grid, startNode=(0, 0), goalNode=(2, 2),
            >>>               heuristic='manhattan', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=None)
            >>> # returning the way back
            >>> astar.getWayBack(astar.mazeStr[2, 2])
            array([[2, 1, 0, 0],
                   [2, 1, 1, 0]])
            >>> astar.getWayBack(astar.mazeStr[1, 2])
            array([[1, 0, 0],
                   [2, 1, 0]])
            >>> astar.getWayBack(astar.mazeStr[1, 1])
            ValueError: Input node is a block. It doesn't have a wayback
        '''
        import numpy as np

        if node.val == -1:
            raise ValueError("Input node is outside the allowed cells")
        elif node.val == 1:
            raise ValueError(
                    "Input node is a block. It doesn't have a wayback")
        pathX = list()
        pathY = list()
        if node.father is not None:
            while node.father != node.pos:
                pathX.append(node.pos[1])
                pathY.append(node.pos[0])
                node = self.mazeStr[node.father]
            pathX.append(self.startNode[1])
            pathY.append(self.startNode[0])
        else:
            pathX = []
            pathY = []
        wayBack = np.array([pathY, pathX])
        return wayBack

    def getPath(self):
        '''Method for obtaining the optimum path between two points into a
        grid-graph through the :math:`\\mathrm{A}^\\ast` algorithm
        `Hart, Nilsson & Raphael (1968) <https://doi.org/10.1109/TSSC.1968.300136>`_.

        Returns:
            (`numpy.ndarray`): :math:`\\left(2 \\times n \\right)` array\
                where :math:`n` is the number of nodes where the path has\
                crossed; the first row of the array contains the abscisses and\
                the second one contains the ordinates of the nodes into the\
                grid-graph

        Examples:
            >>> from numpy import array
            >>> from pybimstab.astar import PreferredPath, Astar
            >>> grid = array([[ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0],
            >>>               [ 0,  0,  0,  0,  1,  1,  0,  0,  0,  0],
            >>>               [ 0,  0,  0,  1,  1,  0,  1,  0,  1,  0],
            >>>               [ 0,  1,  1,  0,  0,  0,  0,  1,  0,  0],
            >>>               [ 0,  0,  0,  1,  0,  0,  0,  1,  1,  1],
            >>>               [ 0,  0,  1,  0,  1,  1,  0,  0,  0,  0],
            >>>               [ 0,  0,  1,  0,  1,  0,  1,  0,  0,  0],
            >>>               [ 1,  1,  0,  1,  0,  0,  0,  0,  1,  1],
            >>>               [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
            >>>               [ 0,  0,  0,  0,  1,  0,  1,  0,  0,  1],
            >>>               [ 0,  1,  0,  0,  0,  1,  0,  0,  0,  0],
            >>>               [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
            >>>               [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
            >>>               [ 1,  0,  0,  0,  1,  0,  0,  0,  0,  0]
            >>>               ])
            >>> coordsIdx = array(
            >>>     [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 13,
            >>>       13, 13, 13, 13, 13, 13, 13],
            >>>      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4,
            >>>       5, 6, 7, 8, 9]])
            >>> preferredPath = PreferredPath(coordsIdx, factor=1)

            >>> # without a forced path
            >>> astar = Astar(grid, startNode=(0, 0), goalNode=(13, 9),
            >>>               heuristic='manhattan', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=None)
            >>> astar.getPath()
            array(
                [[13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  4,  3,  2,  1,  0],
                 [ 9,  9,  9,  9,  8,  8,  7,  7,  6,  5,  4,  3,  2,  1,  0]])

            >>> # with a forced path
            >>> astar = Astar(grid, startNode=(0, 0), goalNode=(13, 9),
            >>>               heuristic='manhattan', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=preferredPath)
            >>> astar.getPath()
            array([
                [13, 13, 13, 13, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                [ 9, 8, 7, 6, 5, 4,  3, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0]])

            >>> from numpy import array
            >>> from pybimstab.astar import PreferredPath, Astar
            >>> grid = array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
                              [1, 1, 0, 0, 1, 0, 1, 0, 0, 0],
                              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
                              [0, 0, 1, 1, 0, 1, 1, 1, 0, 0],
                              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
                              [0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                              [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
                              [0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])
            >>> astar = Astar(grid, startNode=(0, 0), goalNode=(9, 9),
            >>>               heuristic='manhattan', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=None)
            >>> astar.getPath()
            array([[9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                   [9, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0]])
            >>> astar = Astar(grid, startNode=(0, 0), goalNode=(9, 9),
            >>>               heuristic='euclidean', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=None)
            >>> astar.getPath()
            array([[9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0],
                   [9, 8, 7, 6, 5, 4, 5, 5, 4, 3, 2, 1, 0]])
        '''
        import numpy as np

        # Check possible errors in the goal and start cells.
        self.__checkerrors()
        # Required lists
        openList = list()  # open nodes for potential paths
        closeList = list()  # open nodes that have been ruled out
        openList.append(self.mazeStr[self.startNode])
        continueWhile = True
        # Main loop
        while continueWhile:
            # If there is not a solution, stop the program.
            if len(openList) == 0:
                raise ValueError('There is not a path between the start ' +
                                 'and goal nodes')
            # Get the cost of openList
            costList = [node.getGcost(fatherNode=self.mazeStr[node.father]) +
                        node.getHcost(goalNode=self.mazeStr[self.goalNode],
                                      heuristic=self.heuristic,
                                      preferredPath=self.preferredPath)
                        for node in openList]
            # Extract tail-node of current optimum path (currentNode)
            minCostIdx = costList.index(min(costList))  # Index in the list
            currentNode = openList[minCostIdx]  # index of the current node
            openList.pop(minCostIdx)  # delete currentNode from openList
            closeList.append(currentNode)  # add currentNode to closeList
            neighbours = self.getNeighbours(currentNode)  # get the neighbours
            # working on the neighbours of the current node
            for neighbour in neighbours:
                neighbourStr = self.mazeStr[neighbour]
                if (neighbourStr.gCost is not None and
                        np.isinf(neighbourStr.gCost)) or \
                        neighbourStr in closeList:
                    continue
                # verifying if some neighbour is the goal point
                if neighbour == self.goalNode:
                    self.mazeStr[self.goalNode].father = currentNode.pos
                    self.mazeStr[neighbour].getGcost(
                        fatherNode=self.mazeStr[currentNode.pos])
                    self.mazeStr[neighbour].getHcost(
                        goalNode=self.mazeStr[self.goalNode],
                        heuristic=self.heuristic,
                        preferredPath=self.preferredPath)
                    continueWhile = False
                    break
                # If any neighbour is the goal node, A Star continues.
                else:
                    # If the neighbour is in openList: if true, compare the
                    # old gCost (from the previous father) and new gCost (from
                    # the current point which could be the new father) if new
                    # cost is bigger than the old one, leave the original
                    # father, else, change it to the new one.
                    if neighbourStr in openList:
                        oldGcost = neighbourStr.gCost
                        oldFather = neighbourStr.father
                        # changing the attributes: father is current node.
                        self.mazeStr[neighbour].father = currentNode.pos
                        newGcost = self.mazeStr[neighbour].getGcost(
                                fatherNode=self.mazeStr[currentNode.pos])
                        if newGcost > oldGcost:
                            # Coming back to the previous attributes.
                            self.mazeStr[neighbour].father = oldFather
                            self.mazeStr[neighbour].getGcost(
                                    fatherNode=self.mazeStr[oldFather])
                    # If the neighbour is not in the closeList, put it in the
                    # openList.
                    elif neighbourStr not in closeList:
                        self.mazeStr[neighbour].father = currentNode.pos
                        self.mazeStr[neighbour].getGcost(
                            fatherNode=self.mazeStr[currentNode.pos])
                        self.mazeStr[neighbour].getHcost(
                                goalNode=self.mazeStr[self.goalNode],
                                heuristic=self.heuristic,
                                preferredPath=self.preferredPath)
                        openList.append(self.mazeStr[neighbour])
        optimumPath = self.getWayBack(self.mazeStr[self.goalNode])
        setattr(self, 'optimumPath', optimumPath)
        return optimumPath

    def plot(self, plotPreferredPath=False):
        '''Method for generating a graphic of the optimum path returned from
        the ``astar`` method.

        Args:
            plotForcedPath (`bool`): logical variable to check if the forced
                path exists and is wanted to plot.

        Returns:
            (`matplotlib.figure.Figure`): object with the matplotlib structure\
                of the plot. You might use it to save the figure for example.

        Examples:
            >>> from numpy import array
            >>> from pybimstab.astar import PreferredPath, Astar
            >>> grid = array([[0,  0,  0,  0,  1,  0,  0,  0,  0,  0],
            >>>               [0,  0,  0,  0,  1,  1,  0,  0,  0,  0],
            >>>               [0,  0,  0,  1,  1,  0,  1,  0,  1,  0],
            >>>               [0,  1,  1,  0,  0,  0,  0,  1,  0,  0],
            >>>               [0,  0,  0,  1,  0,  0,  0,  1,  1,  1],
            >>>               [0,  0,  1,  0,  1,  1,  0,  0,  0,  0],
            >>>               [0,  0,  1,  0,  1,  0,  1,  0,  0,  0],
            >>>               [1,  1,  0,  1,  0,  0,  0,  0,  1,  1],
            >>>               [1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
            >>>               [0,  0,  0,  0,  1,  0,  1,  0,  0,  1],
            >>>               [0,  1,  0,  0,  0,  1,  0,  0,  0,  0],
            >>>               [1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
            >>>               [0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
            >>>               [1,  0,  0,  0,  1,  0,  0,  0,  0,  0]
            >>>               ])
            >>> coordsIdx = array(
            >>>     [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 13,
            >>>       13, 13, 13, 13, 13, 13, 13],
            >>>      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4,
            >>>       5, 6, 7, 8, 9]])
            >>> preferredPath = PreferredPath(coordsIdx, factor=1)
            >>> for typePP in [None, preferredPath]:
            >>>     astar = Astar(grid, startNode=(0, 0), goalNode=(13, 9),
            >>>                   heuristic='manhattan', reverseLeft=True,
            >>>                   reverseUp=True, preferredPath=typePP)
            >>>     fig = astar.plot(plotPreferredPath=True)

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/astar_example1a.svg
                :alt: astar_example1b

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/astar_example1b.svg
                :alt: astar_example1a

            .. only:: html

               :download:`example script<../examples/figuresScripts/astar_example1.py>`.

            >>> from numpy import array
            >>> from pybimstab.astar import Astar
            >>> grid = array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
                              [1, 1, 0, 0, 1, 0, 1, 0, 0, 0],
                              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
                              [0, 0, 1, 1, 0, 1, 1, 1, 0, 0],
                              [0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
                              [0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                              [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
                              [0, 0, 1, 0, 0, 1, 0, 0, 0, 0]])
            >>> for heuristic in ['manhattan', 'euclidean']:
            >>>     astar = Astar(grid, startNode=(0, 0), goalNode=(9, 9),
            >>>                   heuristic=heuristic, reverseLeft=True,
            >>>                   reverseUp=True, preferredPath=None)
            >>>     fig = astar.plot()

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/astar_example2a.svg
                :alt: astar_example1b

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/astar_example2b.svg
                :alt: astar_example2a

            .. only:: html

               :download:`example script<../examples/figuresScripts/astar_example2.py>`.

            >>> from numpy import array
            >>> from pybimstab.bim import BlocksInMatrix
            >>> from pybimstab.astar import Astar
            >>> seed = 111  # for repeatibilty
            >>> boundary = array([[-5, 0, 5, 0, -5], [0, 10, 0, -10, 0]])
            >>> bim = BlocksInMatrix(slopeCoords=boundary, blockProp=0.2,
            >>>                      tileSize=0.4, seed=seed)
            >>> astar = Astar(bim.grid, startNode=(0, 12), goalNode=(49, 12),
            >>>               heuristic='manhattan', reverseLeft=True,
            >>>               reverseUp=True, preferredPath=None)
            >>> fig = astar.plot()

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/astar_example3.svg
                :alt: astar_example3

            .. only:: html

               :download:`example script<../examples/figuresScripts/astar_example3.py>`.

        '''
        from matplotlib import pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap as newcmap
        import numpy as np

        # Variables to control the color map and its legend
        m, n = self.grid.shape  # dimension of the grid
        yCells, xCells = np.mgrid[slice(-0.5, m-0.5+1, 1),  # grid coordinates
                                  slice(-0.5, n-0.5+1, 1)]
        if np.any(self.grid == -1):  # colormap
            cmap = newcmap.from_list(
                    'BIMcmap', ['white', 'lightgray', 'black'], 3)
            ticks = [-1+0.333, 0, 1-0.333]
            ticksLabels = ['Not allowed cells', 'Allowed cells',
                           'Hindered cells']
        else:
            cmap = newcmap.from_list('BIMcmap', ['lightgray', 'black'], 2)
            ticks = [0.25, 0.75]
            ticksLabels = ['Allowed cells', 'Hindered cells']
        if m > 50 or n > 50:
            edgecolor = 'None'
        else:
            edgecolor = 'k'
        # Plot body
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bar = ax.pcolor(xCells, yCells, self.grid, cmap=cmap,
                        edgecolor=edgecolor)
        ax.plot(self.startNode[1], self.startNode[0], '*r', label='Start node')
        ax.plot(self.goalNode[1], self.goalNode[0], '.r', label='Goal node')
        if plotPreferredPath and self.preferredPath is not None:
            ax.plot(self.preferredPath.coordsIdx[1],
                    self.preferredPath.coordsIdx[0],
                    ':r', lw=1.5, label='Preferred path')
        ax.plot(self.optimumPath[1], self.optimumPath[0], '-r', lw=1.5,
                label='Optimum path')
        # Configuring the colorbar
        bar = plt.colorbar(bar, ax=ax, ticks=ticks, pad=0.01,
                           shrink=0.15, aspect=3)
        bar.ax.set_yticklabels(ticksLabels, fontsize='small')
        # Plot settings
        ax.set_aspect(1)
        ax.legend(fontsize='small', bbox_to_anchor=(1.005, 1), loc=2)
        plt.gca().invert_yaxis()  # invert y-axis acording to de grid notation
        fig.tight_layout()
        return fig


# %%
'''
BSD 2 license.

Copyright (c) 2018, Universidad Nacional de Colombia, Exneyder Andres Montoya
    Araque and Ludger O. Suarez-Burgoa.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
