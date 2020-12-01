## yapCAD polyline and polygon sampling and drawing example

print("example5.py -- yapCAD computational geometry and DXF drawing example")
print("create, sample, and draw simple polylines and polygons")

print("""

In this demo we create simple polyline/polygon figures by regularly
sampling circles, then addign those points to a list.  By repeating
the first point as the last point, we close the figure.

We then use drawable objects to draw the figures, and sample each
polyline to create a series of circles drawn on top.  The circles
increase in size to show the direction of increasing sample parameter.
""")

from yapcad.ezdxf_drawable import *
from yapcad.geom import *

#set up DXF rendering
d=ezdxfDraw()
filename="example5-out"
print("Output file name is {}.dxf".format(filename))
d.filename = filename

d.layer = 'DOCUMENTATION'

d.draw_text("yapCAD", point(5,15),\
            attr={'style': 'OpenSans-Bold',
                  'height': 1.5})

d.draw_text("example5.py",
            point(5,12))
d.draw_text("constructing simple polys",
            point(5,10))
d.draw_text("and sampling outlines",
            point(5,8.5))
d.layer = False
# make circles centered at -10.0,10.0
a = point(-10.0,10.0)

# make circles centered at 10.0,-10.0
b = point(10.0,-10.0)

## sample points from the circles to make regular polygons

## this will be a list of polys, which is to say a list of lists of
## points
polys = []

## NOTE: a simple polyline is just a list of points
for i in range(4,7):
    poly = []
    circ1 = arc(a,i**1.4)
    for j in range(i):
        x=j/i
        poly.append(samplearc(circ1,x)) # sample a point and add it
    poly.append(poly[0]) # make the figure closed
    polys.append(poly) # add the finished poly to the list

for i in range(7,10):
    poly = []
    circ2 = arc(b,(i*0.7)**1.4)
    for j in range(i):
        x=j/i
        poly.append(samplearc(circ2,x)) # sample a point and add it
    poly.append(poly[0]) # make the figure closed
    polys.append(poly) # add the finished poly to the list


## draw the geometry list of polys
d.draw(polys)

# Sample spome points along each poly and use the centers to draw
# circles of increasing size.  This allows the visual confirmation of
# even spacing and shows the counter-clockwise progression of samples

d.layer = 'DRILLS' #select drills layer for circles

i = 4
for ply in polys:
    for j in range(i*3):
        p = samplepoly(ply,j/(i*3))
        d.draw(arc(p,0.2+j/(i*10)))
    i = i+1
    
d.display()
