# -*- coding: utf-8 -*-
'''
Module for defining the repetitive functions to manipulate the
`shapely <https://pypi.python.org/pypi/Shapely>`_ objects that are used in the
classes for performing the slope stabity analysis.
'''


# %%
import numpy as np
from shapely.geometry import LineString, Point


def getIntersect(x, y1, y2):
    '''Intersects two lines that have the same abscise points

    Args:
        x (`list` or `tuple`): abscises of both lines
        y1 (`list` or `tuple`): ordinates of first line
        y2 (`list` or `tuple`): ordinates of second line

    Returns:
        (`tuple`): Coordinates of the interesection point

    Examples:
        >>> from pybimstab.tools import getIntersect
        >>> getIntersect(x=[0, 1], y1=[0, 0], y2=[1, -1])
        (0.5, 0.0)
    '''
    l1 = LineString(np.array([x, y1]).T)
    l2 = LineString(np.array([x, y2]).T)
    try:
        pt = l1.intersection(l2)
        if pt.type != 'Point':
            return None
        else:
            return pt.x, pt.y
    except Exception:
        return None

def cutLine(line, distance):
    '''Cuts a line in two parts at a distance measured from its starting point.

    Args:
        line (`LineString`): Any line or polyline from the ``shapely`` module.
        distance (`int` or `float`): The absolute distance measured from the
            starting point of the line.

    Returns:
        (`tuple`): Two lines that merged at the common point conform the\
            original input line.

    Examples:
        >>> from shapely.geometry import LineString
        >>> from pybimstab.tools import cutLine
        >>> line = LineString([(0, 0), (5, 0)])
        >>> line1, line2 = cutLine(line, 2.5)
        >>> list(line1.coords)
        [(0.0, 0.0), (2.5, 0.0)]
        >>> list(line2.coords)
        [(2.5, 0.0), (5.0, 0.0)]
    '''
    coords = list(line.coords)
    for i, pt in enumerate(coords):
        distAtPt = line.project(Point(pt))
        if distAtPt == distance:
            return LineString(coords[:i+1]), LineString(coords[i:])
        if distAtPt > distance:
            interpPt = line.interpolate(distance)
            line1 = LineString(coords[:i] + [(interpPt.x, interpPt.y)])
            line2 = LineString([(interpPt.x, interpPt.y)] + coords[i:])
            return line1, line2


def upperLowerSegm(line1, line2):
    '''Splits two polylines ath their intersection points and join the upper
    segments in a new polyline and the lower segments in other

    Args:
        line1 (`LineString`): Any line or polyline from the ``shapely`` module.
        line2 (`LineString`): Any line or polyline from the ``shapely`` module.

    Returns:
        (`tuple`): Two lines where each one is given in a (n, 2) array with \
            the coordinates. The first one is the outcome of joinning the\
            upper segments and the second is the outcome of joinning the\
            lower segments.

    Examples:
        >>> from shapely.geometry import LineString
        >>> from pybimstab.tools import upperLowerSegm
        >>> line1 = LineString([(0, 0), (5, 0)])
        >>> line2 = LineString([(0, 5), (5, -5)])
        >>> upperLine, lowerLine = upperLowerSegm(line1, line2)
        >>> upperLine, lowerLine
        (array([[ 0. ,  0. ],
                [ 2.5,  0. ],
                [ 2.5,  0. ],
                [ 5. , -5. ]]),
         array([[0. , 5. ],
                [2.5, 0. ],
                [2.5, 0. ],
                [5. , 0. ]]))
    '''
    # Verification of possible situations
    points = line1.intersection(line2)
    correctLines = True
    # Multiple intersection points (IP)
    if points.geom_type == 'MultiPoint':
        # Two IP; the second is the end of the line2
        if len(points) == 2 and \
                abs(line2.project(points[-1]) - line2.length) < 1e-2:
            points = points[:-1]
        # > 2 IP; and the the end of the line2 is above the line 1
        elif len(points) % 2:
            points = points[:-1]
    # Just one IP and the end of line2 is far enough from line1
    elif points.geom_type == 'Point' and \
            line1.distance(line2.interpolate(line2.length)) > 1e-2:
        points = [points]
    # There is not any real IP; the end of line2 is close enough to line1
    elif line1.distance(line2.interpolate(line2.length)) < 1e-2:
        correctLines = False
        wtEndDist = line2.interpolate(line2.length)
        terrainSlipted = cutLine(line1, line1.project(wtEndDist))
        wt2complete = np.array(terrainSlipted[1].coords)
        waterTabCoords = np.vstack((np.array(line2.coords),
                                    wt2complete))
        upperLine, lowerLine = line1, LineString(waterTabCoords)
    else:
        correctLines = False
        upperLine, lowerLine = line1, line2
    # Make corrections if it is necessary
    if correctLines:
        l1Segmented, l2Segmented = list(), list()
        # Cut the lines in the intersection points
        for pt in points:
            # Manipulation of line 1
            distPt1 = line1.project(pt)  # distance for each interserction
            line1Slipted = cutLine(line1, distPt1)
            line1 = line1Slipted[1]
            l1Segmented.append(np.array(line1Slipted[0]))
            # Manipulation of line 2
            distPt2 = line2.project(pt)  # distance for each interserction
            line2Slipted = cutLine(line2, distPt2)
            line2 = line2Slipted[1]
            l2Segmented.append(np.array(line2Slipted[0]))
        l1Segmented.append(np.array(line1Slipted[1]))
        l2Segmented.append(np.array(line2Slipted[1]))
        # Join the respective segments
        upperLine, lowerLine = l1Segmented[0], l2Segmented[0]
        for segment in range(1, len(points)+1):
            if segment % 2 == 1:
                upperLine = np.vstack((upperLine, l2Segmented[segment]))
                lowerLine = np.vstack((lowerLine, l1Segmented[segment]))
            else:
                upperLine = np.vstack((upperLine, l1Segmented[segment]))
                lowerLine = np.vstack((lowerLine, l2Segmented[segment]))
    return upperLine, lowerLine


def getPointAtX(line, x):
    '''Interesects a vertical line at the abscice ``x`` with the input line and
    retuns the intersection point.

    The input line must not have more than one intersection with the vertical
    line.

    If the abscise is out of the horizontal range of the line, then the nearest
    end is returned.

    Args:
        line (`LineString`): Any line or polyline from the ``shapely`` module.
        x (`int` or `float`): The abscise where the vertical line is wanted.

    Returns:
        (`tuple`): Two lines where each one is given in a (n, 2) array with \
            the coordinates. The first one is the outcome of joinning the\
            upper segments and the second is the outcome of joinning the\
            lower segments.

    Examples:
        >>> from shapely.geometry import LineString
        >>> from pybimstab.tools import upperLowerSegm, getPointAtX
        >>> line = LineString([(0, 0), (5, 0)])
        >>> x = 2.5
        >>> point = getPointAtX(line, x)
        >>> point.x, point.y
        (2.5, 0.0)
    '''
    xMin, yMin, xMax, yMax = np.array(line.bounds)
    vertLine = LineString([(x, yMin - 1), (x, yMax + 1)])
    if xMin > x:
        point = line.interpolate(0, normalized=True)
    elif x > xMax:
        point = line.interpolate(1, normalized=True)
    else:
        point = line.intersection(vertLine)
    if point.geom_type is not 'Point':
        point = point.centroid
    return point


def extractSegment(line, distance1, distance2):
    '''Extracts a segment from a polyline (shapely LineString) given two
    distances measured from its starting point.

    The following rule must be satisfied: `distance1 > distance2`

    Args:
        line (`LineString`): Any line or polyline from the ``shapely`` module.
        distance1 (`int` or `float`): The first distance measured from the
            starting pint of the line.
        distance2 (`int` or `float`): The second distance measured from the
            starting pint of the line.

    Returns:
        (`tuple`): Two lines where each one is given in a (n, 2) array with \
            the coordinates. The first one is the outcome of joinning the\
            upper segments and the second is the outcome of joinning the\
            lower segments.

    Examples:
        >>> from shapely.geometry import LineString
        >>> from pybimstab.tools import upperLowerSegm, extractSegment
        >>> line = LineString([(0, 0), (5, 0)])
        >>> distance1, distance2 = 1, 3
        >>> segment = extractSegment(line, distance1, distance2)
        >>> list(segment.coords)
        [(1.0, 0.0), (3.0, 0.0)]
    '''
    # Split the original line at the first distance
    if distance1 == 0:
        rightSegment = line
    else:
        __, rightSegment = cutLine(line, distance1)
    if distance2 >= line.length:
        segment = rightSegment
    else:  # Split the line2 at the length equal to distance2 - distance1
        segment, __ = cutLine(rightSegment, distance2 - distance1)
    return segment


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
