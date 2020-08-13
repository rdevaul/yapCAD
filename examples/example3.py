## drawing and computational geometry example for yapCAD
print("example3.py -- yapCAD computational geometry and DXF drawing example")
print('''
Create three circles, then compute intersections between circles and
draw tangent lines between pairs of circles.''')

from ezdxf_drawable import *
from geom import *

#set up DXF rendering
d=ezdxfDraw()

filename="example3-out"
print("\nOutput file name is {}.dxf".format(filename))
d.saveas(filename)

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

d.draw(circ1)
d.draw(circ2)
d.draw(circ3)

d.draw(tl1[0])
d.draw(tl1[1])

d.draw(tl2[0])
d.draw(tl2[1])

d.pointstyle='o'
for p in points:
    d.draw(p)

d.display()
