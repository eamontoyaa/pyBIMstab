from shapely.geometry import LineString
from pybimstab.slices import MaterialParameters, SliceStr
material = MaterialParameters(cohesion=15, frictAngle=23,
                              unitWeight=17)
terrainLS = LineString([(6, 8), (7, 6.5)])
slipSurfLS = LineString([(6, 3.395), (7, 2.837)])
watertabLS = LineString([(6, 5), (7, 4)])
slice_ = SliceStr(material, terrainLS, slipSurfLS, watertabLS,
                  bim=None)
fig = slice_.plot()
