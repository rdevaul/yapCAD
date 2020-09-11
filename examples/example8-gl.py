## yapCAD poly() intersecton and drawing examples

print("example8-gl.py -- yapCAD computational geometry openGL drawing example")
print("""

This is a demo of the Polygon() and mirror capabilities of yapCAD,
making and then mirroring random "flowers."
""")

from geom import *
from poly import *
import random

from example8 import *

if __name__ == "__main__":

    from pyglet_drawable import *
    dd=pygletDraw()
    dd.magnify = 0.01
    glist = mirrorArray()

    #drawLegend(dd)
    dd.set_linecolor('white')
    dd.draw(glist)
    dd.display()
