## drawing and computational geometry example for yapCAD
print("example3.py -- yapCAD computational geometry and DXF drawing example")
print('''
Create three circles, then compute intersections between circles and
draw tangent lines between pairs of circles.''')

from yapcad.ezdxf_drawable import *
from yapcad.geom import *

#set up DXF rendering
d=ezdxfDraw()

filename="example3-out"
print("\nOutput file name is {}.dxf".format(filename))
d.filename = filename

d.layer = 'DOCUMENTATION' # select the DOCUMENTATION layer
d.draw_text("yapCAD", point(-9,7),\
            attr={'style': 'OpenSans-Bold',
                  'height': 1.5})

d.draw_text("example3.py",
            point(-9,5.0))
d.draw_text("circles, tangents, and intersections",
            point(-9,3.5))
d.layer = False # select the default layer

a = point(5,5)
b = point(-5,5)
c = point(0,0)

circ1 = arc(a,5)
circ2 = arc(b,7.5)
circ3 = arc(c,1)

tl1 = circleCircleTangentsXY(circ1,circ2)
tl2 = circleCircleTangentsXY(circ3,circ1)

i1 = intersectXY(circ1,circ2,False)
i2 = intersectXY(circ2,circ3,False)

points = i1 + i2

d.draw(circ1)
d.draw(circ2)
d.draw(circ3)

d.draw(tl1[0])
d.draw(tl1[1])

d.draw(tl2[0])
d.draw(tl2[1])

## configure point drawing, draw on DOCUMENTATION layer
d.pointstyle='o'
d.polystyle='points'
d.layer = 'DOCUMENTATION'
d.draw(points)

d.display()
