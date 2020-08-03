## First drawing example for yapCAD

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
drawable=ezdxfDraw()

## make dxf-renderable geometry

# make a point located at 10,10 in the x-y plane, rendered as a small
# cross and circle
point=Point(vect(10,10),"xo")

# make a line segment between the points -5,10 and 10,-5 in the x-y plane
line=Line(vect(-5,10),
          vect(10,-5))

# make an arc with a center at 0,3 with a radius of 3, from 45 degrees
# to 135 degrees
arc=Arc(vect(0,3),3,45,135)


point.draw()
line.draw()
arc.draw()

drawable.saveas("example1-out")
