## drawing and computational geometry example for yapCAD
print("example3.py -- yapCAD computational geometry and DXF drawing example")
print('''
Create three circles, then compute intersections between circles and
draw tangent lines between pairs of circles.''')

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
drawable=ezdxfDraw()

filename="example3-out"
print("\nOutput file name is {}.dxf".format(filename))
drawable.saveas(filename)

a = point(5,5)
b = point(-5,5)
c = point(0,0)

circ1 = arc(a,5)
circ2 = arc(b,7.5)
circ3 = arc(c,1)

tl1 = circleCircleTangentsXY(circ1,circ2)
tl2 = circleCircleTangentsXY(circ3,circ1)

i1 = arcArcIntersectXY(circ1,circ2,False)
i2 = arcArcIntersectXY(circ2,circ3,False)

points = i1 + i2

arc1 = Arc(circ1)
arc2 = Arc(circ2)
arc3 = Arc(circ3)

tl11 = Line(tl1[0])
tl12 = Line(tl1[1])

tl21 = Line(tl2[0])
tl22 = Line(tl2[1])

for p in points:
    P = Point(p,'o')
    P.draw()

arc1.draw()
arc2.draw()
arc3.draw()

tl11.draw()
tl12.draw()

tl21.draw()
tl22.draw()

drawable.display()
