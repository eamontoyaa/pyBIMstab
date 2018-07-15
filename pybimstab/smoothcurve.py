# -*- coding: utf-8 -*-
'''
Module for defining the class related to the curve softener.
'''


# %%
class SmoothCurve:
    '''Creates an instance of an object that defines a curve smoother than the
    input through the :math:`k`-order B-Spline method for interpolation. ::

        SmoothCurve(x, y, k=3, n=300)

    Attributes:
        x (`tuple`, `list` or `numpy.ndarray`): abscisses of the curve to
            smooth.
        y (`tuple`, `list` or `numpy.ndarray`): ordinates of the curve to
            smooth. It must have the same length of ``x``.
        k (`int`): interpolation order.
        n (`int`): number of points of the returned smooth curve

    Note:
        The class ``SmoothCurve`` requires `numpy <http://www.numpy.org/>`_,
        `matplotlib <https://matplotlib.org/>`_. and
        `scipy <https://scipy.org/>`_.

    Examples:
        >>> from pybimstab.smoothcurve import SmoothCurve
        >>> x = [9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 1, 0, 0, 0, 0]
        >>> y = [9, 8, 7, 7, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        >>> curve = SmoothCurve(x=x, y=y)
        >>> curve.__dict__.keys()
        dict_keys(['x', 'y', 'k', 'n', 'smoothing'])
        '''

    def __init__(self, x, y, k=3, n=300):
        '''
        SmoothCurve(x, y, k=3, n=300)
        '''
        from numpy import array
        self.x = array(x)
        self.y = array(y)
        self.k = k
        self.n = n
        # Smoothing the curve
        self.smooth()

    def smooth(self):
        '''Method to generate a smooth curve from the points input through the
        :math:`k`-order B-Spline method.

        Returns:
            (`numpy.ndarray`): :math:`\\left(2 \\times n \\right)` array\
                where :math:`n` is the number of nodes where the path has\
                crossed; the first row of the array contains the abscisses and\
                the second one contains the ordinates of the nodes into the\
                grid-graph.

        Examples:
            >>> from pybimstab.smoothcurve import SmoothCurve
            >>> x = [9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 1, 0, 0, 0, 0]
            >>> y = [9, 8, 7, 7, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0]
            >>> curve = SmoothCurve(x=x, y=y, n=10)
            >>> curve.smooth()
            array([[9.        , 7.56984454, 6.11111111, 4.66666667, 3.2222222,
                    1.77960677, 1.00617284, 0.77206219, 0.01463192, 0.      ],
                   [9.        , 7.05749886, 7.77206219, 8.        , 7.9215821,
                    6.77777778, 5.33333333, 3.88888889, 2.43015546, 0.      ]])
        '''
        import numpy as np
        from scipy.interpolate import splev

        # Defining the knots vector, with k ending equal knots.
        length = len(self.x)
        t = np.linspace(0, 1, length-self.k+1, endpoint=True)
        t = np.append(np.zeros(self.k), t)
        t = np.append(t, np.ones(self.k))
        # Sequence of length 3 containing the knots, coefficients, and degree
        # of the spline to pass it as the tck argument to splev, the function
        # that will evaluate the curve.
        tck = [t, [self.x, self.y], self.k]
        # Required array of the values of the parameter.
        u = np.linspace(0, 1, self.n, endpoint=True)
        # evaluation
        smoothing = np.array(splev(u, tck))
        # Setting the attribute to the instanced object.
        setattr(self, 'smoothing', smoothing)
        return smoothing

    def plot(self):
        '''Method for generating a graphic of the ``smooth`` method output.
        It allows visually to compare the no smoothed line and its smooothed
        version.

        Returns:
            (`matplotlib.figure.Figure`): object with the matplotlib structure\
                of the plot. You might use it to save the figure for example.

        Examples:
            >>> from pybimstab.smoothcurve import SmoothCurve
            >>> x = [9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 1, 0, 0, 0, 0]
            >>> y = [9, 8, 7, 7, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0]
            >>> curve = SmoothCurve(x, y, 1)
            >>> fig = curve.plot()

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/smoothcurve_example1.svg
                :alt: smoothcurve_example1

            .. only:: html

                :download:`example script<../examples/figuresScripts/smoothcurve_example1.py>`.

            >>> import numpy as np
            >>> from pybimstab.smoothcurve import SmoothCurve
            >>> x = np.linspace(0, 2*np.pi, 50)
            >>> y = np.sin(x) + np.random.random(50) * 0.5
            >>> for k in [0, 2, 15]:
            >>>     curve = SmoothCurve(x, y, k)
            >>>     fig = curve.plot()

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/smoothcurve_example2a.svg
                :alt: smoothcurve_example2a

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/smoothcurve_example2b.svg
                :alt: smoothcurve_example2b

            .. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/smoothcurve_example2c.svg
                :alt: smoothcurve_example2c

            .. only:: html

                :download:`example script<../examples/figuresScripts/smoothcurve_example2.py>`.
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # plot original points
        ax.plot(self.x, self.y, 'k--', lw=0.5, marker='x',
                label='Original line')
        # plot smoothed line
        ax.plot(self.smoothing[0], self.smoothing[1], 'k', lw=1.5,
                label='Smoothed curve')
        ax.grid(True, ls='--', lw=0.5)
        ax.legend()
        ax.axis('equal')
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
