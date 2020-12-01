## yapCAD Polyline and Polygon example

print("example6.py -- yapCAD computational geometry and DXF drawing example")
print("create, sample, intersect, and draw Polygons, lines and circles")

print("""
In this demo we sample circles to find the points to create regular
polygons, then we use those points to create circular corners for
Polygon instances.  We draw those polygons in the 'PATHS' layer.

We also demonstrate calculating the intersection of a line with the
polygons in both parameter space and in XY space, and draw the
intersection points in the 'DOCUMENTATION' layer.

Finally, we sample each Polygon instance to create a series of
evenly-spaced circular 'drill holes' along the contour of the poly.
The circles increase in size to show the direction of increasing
sample parameter.  These are rendered in the 'DRILLS' layer.
""")

from yapcad.ezdxf_drawable import *
from yapcad.geom import *
from yapcad.poly import *

#set up DXF rendering
d=ezdxfDraw()
filename="example6-out"
print("Output file name is {}.dxf".format(filename))
d.filename = filename

## Put some documentary text on the drawing
d.layer = 'DOCUMENTATION'

d.draw_text("yapCAD", point(5,15),\
            attr={'style': 'OpenSans-Bold',
                  'height': 1.5})

d.draw_text("example6.py",
            point(5,12))
d.draw_text("polygons, circles, lines,",
            point(5,10))
d.draw_text("and intersections",
            point(5,8.5))
d.layer = False # back to default layer

# make circles centered at -10.0,10.0
a = point(-10.0,10.0)

# make circles centered at 10.0,-10.0
b = point(10.0,-10.0)

## We will store instances of the Polygon() class in this list
polys1 = []

## select the 'PATHS' drawing layer
d.layer = 'PATHS'

## create the first group of concentric rounded polys
for i in range(4,7): # number of sides, 4 to 6
    poly = Polygon()
    polys1.append(poly)
    ## this is the circle we will sample to get the corners of the poly
    circ1 = arc(a,i**1.4-1.0)
    for j in range(i):
        x=j/i
        p = samplearc(circ1,x) ## this is the location of the corner
        poly.addArc(arc(p,1.0)) ## add a circle of radius 1 as the corner

## create the second group of concentric rounded polys
for i in range(7,10): # number of sides, 7 to 9
    poly = Polygon()
    polys1.append(poly)
    circ2 = arc(b,(i*0.7)**1.4-1.0)
    for j in range(i):
        x=j/i
        poly.addArc(arc(samplearc(circ2,x),1.0))

## create a line between the centers of the two groups, and draw it
l = line(a,b)
d.draw(l)

## now, calculate the intersections of this line with all of the polys
## we've just created. We will do that in line parameter space as well
## as simply returning the intersection points.
pp = []
uu1s = []
for ply in polys1:
    z = intersectXY(l,ply.geom(),inside=True)
    u = intersectXY(l,ply.geom(),inside=False,params=True)
    if z:
        pp = pp + z
    if u:
        uu1s = uu1s + u[0]
    # also, draw that poly while we are at it
    d.draw(ply.geom())


## for every intersection point on the line itself, draw it as an 'x'
## in the DOCUMENTATION layer

d.polystyle='points'
d.pointstyle='x'
d.layer = 'DOCUMENTATION'
d.draw(pp)

## for every intersection in line parameter space, including
## intersections outside the 0 <= u <= 1 interval, sample that point
## and draw it as an 'o'.  This means that for the intersection points
## inside the line interval, they will be drawn as a superimposed 'o'
## on the previously drawn 'x'.

d.pointstyle='o'
for u in uu1s:
    p=sampleline(l,u)
    d.draw(p)

## Now, let's draw some evenly-spaced "drill holes" of successivly
## increasing size, tracing the outline of our polys.
d.layer = 'DRILLS' 
i = 4
for ply in polys1:
    for j in range(i*3):
        p = ply.sample(j/(i*3))
        d.draw(arc(p,0.2+j/(i*10)))
    i = i+1

d.layer = False ## Reset to the default layer -- not really necessary
d.display() ## Send the drawing to the file
