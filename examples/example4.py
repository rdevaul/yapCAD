## drawing and computational geometry example for yapCAD
print("example4.py -- yapCAD computational geometry and DXF drawing example")
print("create and render rounded polygons")
print('''

This example demonstrates how to work with the yapCAD Polygon class.
The Polygon class is a generalization of the simple multi-vertex
polygon (see the poly() function in geom.py) intended to support the
creation of 'rounded' geometry and to allow the calculation of
iso-distance lines.

In this case, we create a regular pentagon, then grow that by a
distance of 2.0 to create a rounded version that represents an
iso-distance line.  Then grow it again by 4 to crate another poly that
we sample to find the centers for 11 circles that are equally spaced
along the isodistance line d = 4.0''')

from yapcad.ezdxf_drawable import *
from yapcad.geom import *
from yapcad.poly import *

#set up DXF rendering
d=ezdxfDraw()

filename="example4-out"
print("\nOutput file name is {}.dxf".format(filename))
d.filename = filename

## Put some documentary text on the drawing
d.layer = 'DOCUMENTATION'

d.draw_text("yapCAD", point(5,15),\
            attr={'style': 'OpenSans-Bold',
                  'height': 1.5})

d.draw_text("example4.py",
            point(5,12))
d.draw_text("polygons, circles, and",
            point(5,10))
d.draw_text("geometry sampling",
            point(5,8.5))
d.layer = False

# make a circle with radius 10, centered at the origin
a = point(0,0)
circ = arc(a,10.0)

# sample 5 points from the circle
points = []
for i in range(5):
    x=i/5.0
    points.append(samplearc(circ,x))

# make a circle of radius 2 at the center
center = arc(a,2.0)
d.draw(center)

# make a polygon and add five points as corners
poly0 = Polygon()
for p in points:
    poly0.addPoint(p)

# make a polygon and add five circles with radius 2.0 as corners
poly1 = Polygon()
for p in points:
    a = arc(p,2.0)
    poly1.addArc(a)

# make a third polygon and add circles with radius 4.0 as corners
poly2 = Polygon()
for p in points:
    a = arc(p,4.0)
    poly2.addArc(a)

# draw the generated geometry
d.draw(poly0.geom())
d.draw(poly1.geom())

# sample 12 points from a polygon
for i in range(12):
    x = i/12.0
    p = poly2.sample(x)
    # add some points to the drawing
    if i == 0:
        d.pointstyle = 'xo'
    else:
        d.pointstyle = 'o'
    d.layer = 'DOCUMENTATION'
    d.draw(p)
        
    d.layer = 'DRILLS'
    d.draw(arc(p,2.0))

d.display()
