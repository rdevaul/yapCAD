## yapCAD poly() intersecton and drawing examples

print("example8-gl.py -- yapCAD computational geometry openGL drawing example")
print("""
This is a demo of the Polygon() and mirror capabilities of yapCAD,
making and then mirroring random "flowers."
""")

from yapcad.geom import *
from yapcad.poly import *
import random

from example8 import *

if __name__ == "__main__":

    from yapcad.pyglet_drawable import *
    dd=pygletDraw()
    dd.cameradist=150.0
    glist = mirrorArray()

    #drawLegend(dd)
    dd.set_linecolor('white')
    dd.draw(glist)
    dd.display()
