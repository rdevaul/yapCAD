## yapCAD poly() intersecton and drawing examples

print("example7.py -- yapCAD computational geometry and DXF drawing example")
print("""
In this demo we create open and closed polyline figures and calculate
intersections with lines and arcs.  We draw the elememnts so we can
visibly confirm the correctness of the intersection calculations.""")

from yapcad.ezdxf_drawable import *
from yapcad.geom import *

#set up DXF rendering
dd=ezdxfDraw()
filename="example7-out"
print("Output file name is {}.dxf".format(filename))
dd.filename = filename

dd.layer = 'DOCUMENTATION'

dd.draw_text("yapCAD", point(-12,13),\
            attr={'style': 'OpenSans-Bold',
                  'height': 1.5})

dd.draw_text("example7.py",
            point(-12,10))
dd.draw_text("polygons, polylines, arcs,",
            point(-12,8))
dd.draw_text("and intersections",
            point(-12,6.5))
dd.layer = False

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

dd.layer = 'PATHS'
dd.linecolor = 1

dd.draw(line1)
dd.draw(line2)
dd.draw(line3)
dd.draw(line4)

dd.linecolor = 4
dd.draw(arc1)
dd.draw(arc2)
dd.draw(arc3)
dd.draw(arc4)

dd.linecolor = False

## Do some intersection calculation

int0 = intersectXY(line1,pol1,True)
int1 = intersectXY(line2,pol1,True)
int2 = intersectXY(arc1,pol1,True)
int2u =intersectXY(arc1,pol1,params=True)
int3 = intersectXY(arc2,pol1,True)
int3u = intersectXY(arc2,pol1,params=True)

int10 = intersectXY(line3,pol2,True)
int11 = intersectXY(line4,pol2,True)
int12 = intersectXY(arc3,pol2,True)
int12u =intersectXY(arc3,pol2,params=True)
int13 = intersectXY(arc4,pol2,True)
int13u = intersectXY(arc4,pol2,params=True)

## draw intersection points on documentation layer

dd.layer = 'DOCUMENTATION'
dd.polystyle = 'points'
dd.draw(int0 + int1 + int10 + int11)

dd.pointstyle='o'
i = 6
for p in int2 + int12:
    dd.draw(p)
    # dd.draw(arc(p,0.1*i))
    i = i+1

dd.pointstyle='x'
i = 6
for u in int2u[1]:
    p = samplepoly(pol1,u)
    #print("sample point: ",vstr(p))
    dd.draw(p)
    #dd.draw(arc(p,0.1*i))
    i = i+1

i = 6
for u in int12u[1]:
    p = samplepoly(pol2,u)
    dd.draw(p)
    #dd.draw(arc(p,0.1*i))
    i = i+1
    
dd.pointstyle='o'
i = 6
for p in int3 + int13:
    dd.draw(p)
    #dd.draw(arc(p,0.1*i))
    i = i+1

dd.pointstyle='x'
i = 6
for u in int3u[1]:
    p = samplepoly(pol1,u)
    #print("sample point: ",vstr(p))
    dd.draw(p)
    #dd.draw(arc(p,0.1*i))
    i = i+1

i = 6
for u in int13u[1]:
    p = samplepoly(pol2,u)
    dd.layer = 'DOCUMENTATION'
    dd.draw(p)
    dd.layer = 'DRILLS'
    dd.draw(arc(p,0.1*i))
    i = i+1


dd.display()
