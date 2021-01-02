import pytest
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.pyglet_drawable import *
from yapcad.poly import *
from yapcad.combine import *

## unit tests for yapCAD geom.py

class TestSample:
    """test sampling-related operations"""

    def test_sample(self):
        dd = pygletDraw()
        dd.linecolor='silver'
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(0,0),3.0,6.66)
        # "close" the arc with a line segment
        gl.append(line(sample(gl,1.0),
                       sample(gl,0.0)))
        dd.draw(gl)
        # turn the geometry list into a poly

        ply = geomlist2poly(gl,minang=10.0,minlen=0.5)
        ply = translate(ply,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply)
        dd.display()

