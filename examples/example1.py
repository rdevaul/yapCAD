## First drawing example for yapCAD
print("example1.py -- yapCAD DXF drawing example")

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
d=ezdxfDraw()

filename="example1-out"
print("\nOutput file name is {}.dxf".format(filename))
d.saveas(filename)

## Put some documentary text on the drawing
d.layerset('DOCUMENTATION') # select the DOCUMENTATION layer

d.draw_text("yapCAD", point(5,15),\
            attr={'style': 'OpenSans-Bold',
                  'height': 2.0})

d.draw_text("example1.py",
            point(5,12))
d.draw_text("drawing primitives",
            point(5,10))
d.layerset() # select the default layer

## make dxf-renderable geometry

# make a point located at 10,8 in the x-y plane, rendered as a small
# red cross and circle

d.pointstyle = 'xo' # set the point rendering style
d.colorset(1) #set color to red (DXF index color)
d.draw(point(10,8))

d.colorset(7) #set color to white
# make a line segment between the points -5,15 and 10,-5 in the x-y plane
d.draw(line(point(-5,15),
            point(10,-5)))

# make an arc with a center at 1,3 with a radius of 6, from 45 degrees
# to 135 degrees
d.draw(arc(point(1,3),6,45,135))

# write out the geometry as example1-out.dxf
d.display()
