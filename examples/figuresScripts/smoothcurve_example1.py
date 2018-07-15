from pybimstab.smoothcurve import SmoothCurve
x = [9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 1, 0, 0, 0, 0]
y = [9, 8, 7, 7, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0]
curve = SmoothCurve(x, y)
fig = curve.plot()
