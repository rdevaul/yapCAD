## drawing and computational geometry example for yapCAD
print("yapCAD computational geometry and DXF drawing example")
print("create and render rounded polygons")
print("\n\
This example demonstrates how 'rounded' polygons can be used to\n\
calculate standoff distance from conventional polygons.  In this\n\
case, we create a regular pentagon, then grow that by a distance \n\
of 2.0 to create a rounded version that represents an iso-distance \n\
line.  Then grow it again by 4 to crate another poly that we \n\
sample to find the centers for 11 circles that are equally spaced\n\
along the isodistance line d = 4.0")

from ezdxf_drawable import *
from geom import *
from poly import *

#set up DXF rendering
drawable=ezdxfDraw()
drawable.saveas("example4-out")

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
Arc(center).draw()

# make a polygon and add five points as corners
poly0 = Poly()
for p in points:
    poly0.addPoint(p)

# make a polygon and add five circles with radius 2.0 as corners
poly1 = Poly()
for p in points:
    a = arc(p,2.0)
    poly1.addArc(a)

# make a third polygon and add circles with radius 4.0 as corners
poly2 = Poly()
for p in points:
    a = arc(p,4.0)
    poly2.addArc(a)

# create the outlines
poly0.makeoutline()
poly1.makeoutline()
poly2.makeoutline()

def drawOutline(ply):    
    # draw the outline
    for e in ply.outline:
        if isline(e):
            l = Line(e)
            l.draw()
        else: #it's an arc
            a = Arc(e)
            a.draw()
            
drawOutline(poly0)
drawOutline(poly1)
#drawOutline(poly2)

# sample 12 points from a polygon
for i in range(12):
    x = i/12.0
    p = poly2.sample(x)
    # add some points to the drawing
    if i == 0:
        Point(p,'xo').draw()    # mark the first point differently
    else:
        Point(p,'o').draw()
        
    Arc(arc(p,2.0)).draw()

drawable.display()
