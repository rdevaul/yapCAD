## First drawing example for yapCAD
print("example1.py -- yapCAD DXF drawing example")

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
d=ezdxfDraw()

filename="example1-out"
print("\nOutput file name is {}.dxf".format(filename))
d.saveas(filename)

## make dxf-renderable geometry

# make a point located at 10,10 in the x-y plane, rendered as a small
# cross and circle

d.pointstyle = 'xo'
d.draw(point(10,10))

# make a line segment between the points -5,10 and 10,-5 in the x-y plane
d.draw(line(point(-5,10),
                   point(10,-5)))

# make an arc with a center at 0,3 with a radius of 3, from 45 degrees
# to 135 degrees
d.draw(arc(point(0,3),3,45,135))

# write out the geometry as example1-out.dxf
d.display()
