import pytest
import os
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.poly import *
from yapcad.combine import *

## unit tests for yapCAD geom.py

# Control flag for visual tests - set via environment variable or directly
VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw
    
class TestSample:
    """test sampling-related operations"""

    @pytest.mark.visual
    def test_sample(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

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

