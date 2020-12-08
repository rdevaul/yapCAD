## yapCAD boolean operations geometry example

from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *

import math
import random

def makeStars(center,rad=5.0,diam=2.0):
    poly1 = Polygon()
    poly2 = Polygon()

    for i in range(8):
        ang = (pi2/8)*i
        p = poly1
        if i%2 ==0:
            p = poly2
        p.addArc(arc(add(center, point(math.cos(ang)*rad,
                                       math.sin(ang)*rad)),diam))
    return poly1,poly2

def geometry():
    p1 = Polygon()
    x = point(-2,2)
    p1.addArc(arc(x,3.0))

    #p1 = makeRoundRect(4,4,0.2,point(-4.0,4.0))

    p2 = makeRoundRect(7,6,1,point(1,-1))
    b = Boolean('union',[p1,p2])

    p4,p5 = makeStars(point(12,-10))
    b2 = Boolean('union',[p4,p5])

    p6,p7 = makeStars(point(12,10))
    b3 = Boolean('intersection',[p6,p7])

    p8,p9 = makeStars(point(-12,-10))
    b4 = Boolean('difference',[p9,p8])

    p10,p11 = makeStars(point(-12,10))

    return [p1.geom(),p2.geom(),b.geom(),b2.geom(),b3.geom(),b4.geom(),p10.geom(),p11.geom()]
    #return [p1.geom(),p2.geom(),b.geom()]

def testAndDraw(dd):
    dd.linecolor = 'red'

    gl = geometry()

    g1 = gl[0]
    g2 = gl[1]
    g3 = gl[2]
    g4 = gl[3]
    g5 = gl[4]
    g6 = gl[5]
    
    dd.linecolor = 'aqua'

    dd.draw(g3)

    g3r = rotate(g3,45)

    dd.linecolor = 'red'
    dd.draw(g3r)

    dd.linecolor = 'white'
    for i in range(len(g4)):
        dd.linecolor = i%7+1
        dd.draw(g4[i])
    
    dd.linecolor = 'silver'
    dd.draw(translate(g4,point(0,0,-1)))
            
    dd.linecolor = 'white'
    dd.draw(g5)

    dd.linecolor = 'yellow'
    dd.draw(g6)

    dd.linecolor = 'aqua'
    dd.draw(gl[6:])

    dd.display()
    
if __name__ == "__main__":
    import sys
    renderOgl = False
    filename="example11-out"
    oglarg= ("pyglet","opengl","OpenGL")
    dxfarg= ("ezdxf","dxf","DXF")
    print("example11.py -- yapCAD boolean operations demonstration")
    if len(sys.argv) > 1:
        if sys.argv[1] in oglarg:
            renderOgl=True
        elif sys.argv[1] in dxfarg:
            renderOgl=False
        else:
            print(" In this example, we create yapCAD Polygon() geometry generators and")
            print(" perform boolean operations on them to create combined figures.")
            print("syntax: $ python3 {} <rendertype> [filename]".format(sys.argv[0]))
            print("    where <rendertype> is one of {} for OpenGL".format(oglarg))
            print("    or one of {} for DXF".format(dxfarg))
            print("    For DXF, you can optionally specify [filename].dxf as the output file")
            quit()
    if len(sys.argv) > 2 and renderOgl==False:
        filename = sys.argv[2]+".dxf"
    dd = []
    if renderOgl:
        print("OpenGL rendering selected")
        from yapcad.pyglet_drawable import *
        dd=pygletDraw()
    else:
        print("DXF rendering selected")
        from yapcad.ezdxf_drawable import *
        #set up DXF rendering
        dd=ezdxfDraw()
        dd.filename = filename
    print("rendering...")
    testAndDraw(dd)
    print("done")
