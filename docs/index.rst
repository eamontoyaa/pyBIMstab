pyBIMstab
=========

``pybimstab`` is an application software in **Python 3** to evaluate the factor
of safety against sliding of slopes made of Blocks-In-Matrix (BIM) materials. 

The assessment is donde by using the limit equilibrium method through the
General Limit Equilibrium (GLE) method of
`Fredlund & Krahn (1977) <https://doi.org/10.1139/t77-045>`_.

The slip surface has a tortuous geometry and is optimally found by using the
math:`\\mathrm{A}^\\ast` algorithm proposed by 
`Hart, Nilsson & Raphael (1968) <https://doi.org/10.1109/TSSC.1968.300136>`_.

The following plots are the final outcome of two different analysis:

**Homogeneus slope**

.. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/coverPlot1.svg
        :alt: Outcome plot example1

**Slope made of BIM material**

.. figure:: https://rawgit.com/eamontoyaa/pybimstab/master/examples/figures/coverPlot2.svg
        :alt: Outcome plot example2


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   modules
   usage
   authors
   license
   history
   references


Links
=====

* `Documentation <https://pybimstab.readthedocs.io>`_
* `PyPI <https://pypi.org/project/pybimstab>`_
* `GitHub <https://github.com/eamontoyaa/pybimstab>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

License and Copyright
=====================
.. |br| raw:: html

   <br />

Copyright (c) 2018, Universidad Nacional de Colombia, Medellín. |br|
Copyright (c) 2018, Exneyder A. Monotoya-Araque and Ludger O. Suarez-Burgoa. |br|
License BSD-2 or higher.
