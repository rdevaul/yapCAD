## yapCAD poly() intersecton and drawing examples

import sys
from yapcad.geom import *
from yapcad.poly import *
import random

from examples.example8 import *

if __name__ == "__main__":

    print("example8-gl.py -- yapCAD computational geometry openGL drawing example")
    print("""
 This is a demo of the Polygon() and mirror capabilities of yapCAD,
 making and then mirroring random "flowers."  Renders output interactively with OpenGL.
    """)

    try:
        from yapcad.pyglet_drawable import *
        dd=pygletDraw()
    except RuntimeError as err:
        print("\nSkipping OpenGL rendering: {}".format(err))
        sys.exit(0)
    except Exception as err:
        print("\nSkipping OpenGL rendering due to unexpected error: {}".format(err))
        sys.exit(1)
    dd.cameradist=150.0
    glist = mirrorArray()

    def mydrawer(gl):
        g1=[]
        g2=[]
        for g in gl:
            if iscircle(g):
                g1.append(g[0])
            else:
                g2.append(g)
        dd.linecolor = 'white'
        dd.draw(g2)
        dd.linecolor = 'aqua'
        dd.polystyle='points'
        dd.draw(g1)

    mydrawer(glist)
    dd.display()
