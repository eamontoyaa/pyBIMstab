# -*- coding: utf-8 -*-
'''
Module for defining the class related to a Polygon and to verify what points
are inside it.
'''


# %%
class Polygon:
    '''
    Creates an instance of an object that defines a two dimension polygon. ::

        Polygon(coordinates)

    Attributes:
        coordinates ((n, 2) `numpy.ndarray`): Coordinates of vertices of the
            polygon.

    Note:
        The class `Polygon` requires `numpy <http://www.numpy.org/>`_ and
        `matplotlib <https://matplotlib.org/>`_.

    Examples:
        >>> from numpy import array
        >>> from pybimstab.polygon import Polygon
        >>> coords = array([[0. , 1. , 1. , 0. ], [0. , 0. , 1.5, 1.5]])
        >>> polygon = Polygon(coordinates=coords)
        >>> polygon.__dict__
        {'coordinates': array([[0. , 1. , 1. , 0. ],
                               [0. , 0. , 1.5, 1.5]])}
        '''

    def __init__(self, coordinates):
        '''
        Polygon(coordinates)
        '''
        self.coordinates = coordinates

    def isinside(self, x, y, meshgrid=False, want2plot=False):
        '''Method to know if the point(s) (x, y) is inside of the instanced
        polygon.

        ``x`` and ``y`` could be iterable structures, even, they could define a
        meshgrid.

        Args:
            x (`int`, `float` or (n, ) `numpy.ndarray`): abscissa of
                the point(s) to check if is/are inside the polygon.
            y (`int`, `float` or (n, ) `numpy.ndarray`): ordinate of
                the point(s) to check if is/are inside the polygon.
            meshgrid (`bool`): variable to creck if ``x`` and ``y`` define
                a grid. The default value is ``False``.
            want2plot (`bool`): variable to creck if a plot is wanted.
                The default value is ``False``.

        Returns:
            `numpy.ndarray`: Array of boolean values where `True` means\
                that the point is inside and `False` to points that is outside\
                of the polygon.

        Note:
            The method ``isinside`` does not work properly for points on the
            boundaries of the polygon.

        Examples:
            >>> from numpy import array
            >>> from pybimstab.polygon import Polygon
            >>> coords = array([[0. , 1. , 1. , 0. ], [0. , 0. , 1.5, 1.5]])
            >>> x, y = 0.5, 2
            >>> polygon = Polygon(coordinates=coords)
            >>> polygon.isinside(x=x, y=y, meshgrid=False, want2plot=True)
            array([True], dtype=bool)

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/polygon_example1.svg
                :alt: polygon_example1

            .. only:: html

                :download:`example script<../examples/figuresScripts/polygon_example1.py>`.

            >>> from numpy import array
            >>> from pybimstab.polygon import Polygon
            >>> coords = array([[0, 1, 1, 0.], [0, 0, 1.5, 1.5]])
            >>> x = array([0.3, 0.5, 0.7, 1.2, 1.0])
            >>> y = array([0.6, 1., 1.4, 2.4, 0])
            >>> polygon = Polygon(coordinates=coords)
            >>> polygon.isinside(x=x, y=y, meshgrid=False, want2plot=True)
            array([True, True, True, False, False], dtype=bool)

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/polygon_example2.svg
                :alt: polygon_example2

            .. only:: html

                :download:`example script<../examples/figuresScripts/polygon_example2.py>`.

            >>> from numpy import array
            >>> from pybimstab.polygon import Polygon
            >>> coords = array([[0, 1, 1, 0.], [0, 0, 1.5, 1.5]])
            >>> x = array([0.3, 0.5, 0.7, 1.2, 1.0])
            >>> y = array([0.6, 1., 1.4, 2.4, 0])
            >>> polygon = Polygon(coordinates=coords)
            >>> polygon.isinside(x=x, y=y, meshgrid=True, want2plot=True)
            array([[ True,  True,  True, False,  True],
                   [ True,  True,  True, False,  True],
                   [ True,  True,  True, False,  True],
                   [False, False, False, False, False],
                   [False, False, False, False, False]], dtype=bool)

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/polygon_example3.svg
                :alt: polygon_example3

            .. only:: html

                :download:`example script<../examples/figuresScripts/polygon_example3.py>`.

            >>> import numpy as np
            >>> from pybimstab.polygon import Polygon
            >>> xC = [np.cos(theta) for theta in np.linspace(0, 2*np.pi, 100)]
            >>> yC = [np.sin(theta) for theta in np.linspace(0, 2*np.pi, 100)]
            >>> coords = np.array([xC, yC])
            >>> np.random.seed(123)
            >>> x = np.random.uniform(-1, 1, 100)
            >>> y = np.random.uniform(-1, 1, 100)
            >>> polygon = Polygon(coordinates=coords)
            >>> polygon.isinside(x=x, y=y, meshgrid=False, want2plot=True);

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/polygon_example4.svg
                :alt: polygon_example4

            .. only:: html

                :download:`example script<../examples/figuresScripts/polygon_example4.py>`.
        '''
        import numpy as np
        from matplotlib import pyplot as plt
        from matplotlib.path import Path

        # Definition of the polygon as a path
        pathPolygon = Path(self.coordinates.T)
        # Defining the proper (x, y) structure-array
        x = np.array(x)
        y = np.array(y)
        if meshgrid:
            xGrid, yGrid = np.meshgrid(x, y)
            xyArray = np.dstack((xGrid, yGrid)).reshape((-1, 2))
            inside = pathPolygon.contains_points(xyArray)
        else:
            xyArray = np.vstack((x, y)).transpose()
            inside = pathPolygon.contains_points(xyArray)

        # Plotting
        if want2plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # plot inside points
            ax.plot(xyArray[inside][:, 0], xyArray[inside][:, 1],
                    '.k', label='Inside')
            # plot outside points
            ax.plot(xyArray[np.invert(inside)][:, 0],
                    xyArray[np.invert(inside)][:, 1],
                    'xk', label='Outside')
            # plot the close polygon
            ax.plot(np.hstack((self.coordinates[0], self.coordinates[0, 0])),
                    np.hstack((self.coordinates[1], self.coordinates[1, 0])),
                    'k', label='Polygon')
            ax.grid(True, ls='--', lw=0.5)
            ax.legend()
            ax.axis('equal')
            fig.tight_layout()
        if meshgrid:
            inside = inside.reshape(xGrid.shape)
        return inside


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
