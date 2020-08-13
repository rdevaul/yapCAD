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

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
drawable=ezdxfDraw()
filename="example5-out"
print("Output file name is {}.dxf".format(filename))
drawable.saveas(filename)

# make circles centered at -10.0,10.0
a = point(-10.0,10.0)

# make circles centered at 10.0,-10.0
b = point(10.0,-10.0)

# sample points from the circles to make regular polygons

points1 = []
points2 = []

for i in range(4,7):
    p = []
    circ1 = arc(a,i**1.4)
    for j in range(i):
        x=j/i
        p.append(samplearc(circ1,x))
    p.append(p[0])
    points1.append(p)

for i in range(7,10):
    p = []
    circ2 = arc(b,(i*0.7)**1.4)
    for j in range(i):
        x=j/i
        p.append(samplearc(circ2,x))
    p.append(p[0])
    points2.append(p)

ply1=[]
ply2=[]

for p in points1:
    ply = poly(p)
    ply1.append(ply)

for p in points2:
    ply = poly(p)
    ply1.append(ply)

for ply in ply1:
    for i in range(1,len(ply)):
        l = line(ply[i-1],ply[i])
        Line(l).draw()


i = 4
for ply in ply1:
    for j in range(i*3):
        p = samplepoly(ply,j/(i*3))
        Arc(arc(p,0.2+j/(i*10))).draw()
    i = i+1
    
drawable.display()
