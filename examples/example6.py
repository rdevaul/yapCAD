## yapCAD Polyline and Polygon example

print("example6.py -- yapCAD computational geometry and DXF drawing example")
print("create, sample, and draw Polylines and Polygons")

print("""

In this demo we largely reproduce example5 using the Polygon class --
we sample circles to find the points to create regular polygons, then we use those points to create circular corners for Polygon instances. 

We then use drawable objects to draw the figures, and sample each
Polygon instance to create a series of circles drawn on top.  The
circles increase in size to show the direction of increasing sample
parameter.  """)

from ezdxf_drawable import *
from geom import *
from poly import *

#set up DXF rendering
d=ezdxfDraw()
filename="example6-out"
print("Output file name is {}.dxf".format(filename))
d.saveas(filename)

# make circles centered at -10.0,10.0
a = point(-10.0,10.0)

# make circles centered at 10.0,-10.0
b = point(10.0,-10.0)

# sample points from the circles to make regular polygons

polys1 = []

for i in range(4,7):
    poly = Polygon()
    polys1.append(poly)
    circ1 = arc(a,i**1.4-1.0)
    for j in range(i):
        x=j/i
        poly.addArc(arc(samplearc(circ1,x),1.0))
    poly.makeoutline()

for i in range(7,10):
    poly = Polygon()
    polys1.append(poly)
    circ2 = arc(b,(i*0.7)**1.4-1.0)
    for j in range(i):
        x=j/i
        poly.addArc(arc(samplearc(circ2,x),1.0))
    poly.makeoutline()

    
def drawOutline(ply):    
    # draw the outline
    for e in ply.outline:
        d.draw(e)

for ply in polys1:
    drawOutline(ply)


i = 4
for ply in polys1:
    for j in range(i*3):
        p = ply.sample(j/(i*3))
        d.draw(arc(p,0.2+j/(i*10)))
    i = i+1
    
d.display()
