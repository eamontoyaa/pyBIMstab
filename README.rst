=========
pyBIMstab
=========

.. image:: https://img.shields.io/badge/Made%20with-Python-brightgreen.svg
        :target: https://www.python.org/
        :alt: made-with-python

.. image:: https://img.shields.io/pypi/v/pybimstab.svg
        :target: https://pypi.python.org/pypi/pybimstab
        :alt: PyPI

.. image:: https://img.shields.io/badge/License-BSD%202--Clause-brightgreen.svg
        :target: https://github.com/eamontoyaa/pybimstab/blob/master/LICENSE
        :alt: License

.. image:: https://readthedocs.org/projects/pybimstab/badge/?version=latest
        :target: https://pybimstab.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status



``pybimstab`` is an application software in **Python 3** to evaluate the factor
of safety against sliding of slopes made of Blocks-In-Matrix (BIM) materials. 

The assessment is donde by using the limit equilibrium method through the
General Limit Equilibrium (GLE) method of
`Fredlund & Krahn (1977) <https://doi.org/10.1139/t77-045>`_.

The slip surface has a tortuous geometry and is optimally found by using the
A-star algorithm proposed by 
`Hart, Nilsson & Raphael (1968) <https://doi.org/10.1109/TSSC.1968.300136>`_.

The following plots are the final outcome of two different analysis:

**Homogeneus slope**

.. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/coverPlot1.svg
        :alt: Outcome plot example1

**Slope made of BIM material**

.. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/coverPlot2.svg
        :alt: Outcome plot example2


Features
--------

* `Documentation <https://pybimstab.readthedocs.io>`_
* `PyPI <https://pypi.org/project/pybimstab>`_
* `GitHub <https://github.com/eamontoyaa/pybimstab>`_
* Open source and free software: `BSD-2-Clause <https://opensource.org/licenses/BSD-2-Clause>`_.


Requirements
------------

The code was written in Python 3. The packages `numpy <http://www.numpy.org/>`_,
`scipy <https://www.scipy.org/>`_, `matplotlib <https://matplotlib.org/>`_
and `shapely <https://pypi.org/project/Shapely/>`_ are
required for using ``pybimstab``. All of them are
downloadable from the PyPI repository by opening a terminal and typing the
following code lines:


::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install shapely


Installation
------------


To install ``pybimstab`` open a terminal and type:

::

    pip install pybimstab


Example
-------

To produce the plot shown above execute the following script

::

    from numpy import array
    from pybimstab.slope import NaturalSlope
    from pybimstab.watertable import WaterTable
    from pybimstab.bim import BlocksInMatrix
    from pybimstab.slipsurface import CircularSurface, TortuousSurface
    from pybimstab.slices import MaterialParameters, Slices
    from pybimstab.slopestabl import SlopeStabl
    terrainCoords = array(
        [[-2.49, 0.1, 1.7, 3.89, 5.9, 8.12, 9.87, 13.29, 20.29,
          21.43, 22.28, 23.48, 24.65, 25.17],
         [18.16, 17.88, 17.28, 15.73, 14.31, 13.58, 13, 3.61, 3.61,
          3.32, 2.71, 2.23, 1.21, 0.25]])
    slope = NaturalSlope(terrainCoords)
    bim = BlocksInMatrix(slopeCoords=slope.coords, blockProp=0.25,
                         tileSize=0.4, seed=12345)
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
        slopeCoords=slope.coords, numSlices=15,
        watertabCoords=watertable.coords, bim=bim)
    stabAnalysis = SlopeStabl(slices, seedFS=1, Kh=0, maxLambda=1)
    fig = stabAnalysis.plot()


References
----------
D. G. Fredlund and J. Krahn. Comparison of slope stability methods of analysis.
Canadian Geotechnical Journal, 14(3)(3):429–439, 1977.

P. Hart, N. Nilsson, and B. Raphael. A formal basis for the heuristic
determination of minimum cost path. IEEE Transactions of Systems Science and
Cybernetics, ssc-4(2):100–107, 1968.

