# -*- coding: utf-8 -*-

"""Top-level package for pyBIMstab."""

from .astar import *
from .bim import *
from .polygon import *
from .slices import *
from .slipsurface import *
from .slope import *
from .slopestabl import *
from .smoothcurve import *
from .tools import *
from .watertable import *

__all__ = ['astar', 'bim', 'polygon', 'slices', 'slipsurface', 'slope',
           'slopestabl', 'smoothcurve', 'tools', 'watertable']

__author__ = """Exneyder A. Montoya-Araque & Ludger O. Suarez-Burgoa"""
__email__ = 'eamontoyaa@gmail.com'
__version__ = '0.1.5'
