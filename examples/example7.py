## yapCAD poly() intersecton and drawing examples

print("example7.py -- yapCAD computational geometry and DXF drawing example")
print("""

In this demo we create open and closed polyline figures and calculate
intersections with lines and arcs.  We draw the elememnts so we can
visibly confirm the correctness of the intersection calculations.""")

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
dd=ezdxfDraw()
filename="example7-out"
print("Output file name is {}.dxf".format(filename))
dd.saveas(filename)

# make some points centerd at the origin to be used as vertices for
# our polys
a = point(0,-5)
b = point(5,0)
c = point(0,5)
d = point(-5,0)

## these are offsets that we will for the above points to create two separate
## offset polyline figures
offset1 = point(-4,-4)
offset2 = point(4,4)

## make an open polyline figure in the XY plane from the points we
## defined above, offset by offset1

## Note, I like lambda expressions and map, but there are lots
## of other ways to accomplish this offset.
pol1 = poly(list(map(lambda x: add(x,offset1),[a,b,c,d])))

## make a closed polyline figure in the XY plane from the points we
## defined above, offset by offset2
pol2 = poly(list(map(lambda x: add(x,offset2),[a,b,c,d,a])))
                                               

line1 = line(offset1,add(offset1,point(4,4)))
line2 = line(add(offset1,point(0,-6)),add(offset1,point(0,6)))

arc1 = arc(offset1,4.0,270,90)
arc2 = arc(offset1,4.5,0,360)

line3 = line(offset2,add(offset2,point(4,4)))
line4 = line(add(offset2,point(0,-6)),add(offset2,point(0,6)))

arc3 = arc(offset2,4.0,270,90)
arc4 = arc(offset2,4.5,0,360)

dd.polystyle='lines'
dd.draw(pol1)
dd.draw(pol2)

dd.draw(line1)
dd.draw(line2)
dd.draw(line3)
dd.draw(line4)

dd.draw(arc1)
dd.draw(arc2)
dd.draw(arc3)
dd.draw(arc4)

## Do some intersection calculation

int0 = intersectSimplePolyXY(line1,pol1,True)
int1 = intersectSimplePolyXY(line2,pol1,True)
int2 = intersectSimplePolyXY(arc1,pol1,True)
int2u =intersectSimplePolyXY(arc1,pol1,params=True)
int3 = intersectSimplePolyXY(arc2,pol1,True)

dd.polystyle = 'points'
#dd.draw(int0)
#dd.draw(int1)
dd.draw(int2)
print("intersectSimplePolyXY(arc1,pol1,params=True): ",int2u)
print("intersectSimplePolyXY(arc1,pol1,True): ",vstr(int2))
#dd.draw(int3)

dd.display()
