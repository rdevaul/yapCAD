## yapCAD poly() intersecton and drawing examples

print("example8.py -- yapCAD computational geometry and DXF drawing example")
print("""

This is a demo of the Polygon() and mirror capabilities of yapCAD,
making and then mirroring random "flowers."
""")

from geom import *
from poly import *
import random

def drawLegend(d):
    ## Put some documentary text on the drawing
    d.layerset('DOCUMENTATION')

    d.draw_text("yapCAD", point(5,15),\
                attr={'style': 'OpenSans-Bold',
                      'height': 1.5})

    d.draw_text("example6.py",
                point(5,12))
    d.draw_text("Polygon() flowers, mirrored geometry",
                point(5,10))
    d.layerset() # back to default layer


def flower(center = point(0,0),
           petals = 12,
           minDiam=5.0,
           maxDiam=15,
           minRadius=20,
           maxRadius=40,
           insideRad = 5, returnPoly = False):

    glist= []
    for i in range(petals):
        angle = i*360/petals
        anrad = angle*pi2/360.0

        an2 = (((i+0.5)/petals)%1.0)*360
        an2rad = an2*pi2/360.0
    
        radius = (maxRadius-minRadius)*random.random()+minRadius
        diam = (maxDiam-minDiam)*random.random()+minDiam
        pnt = add(point(cos(anrad)*radius,sin(anrad)*radius),
                  center)
        pnt2 = add(point(cos(an2rad)*insideRad,sin(an2rad)*insideRad),
                   center)
        a = arc(pnt,diam/2)
        glist.append(a)
        glist.append(pnt2)

    p = Polygon(glist)
    #p.makeoutline()
    if returnPoly:
        return p
    else:
        return p.geom()

def mirrorArray(pnt=point(-45,45)):
    glist = flower(pnt)
    bb = bbox(glist)
    ply = [point(bb[0]), point(bb[1][0],bb[0][1]),
           point(bb[1]), point(bb[0][0],bb[1][1])]
    glist = glist + ply
    glist = glist + mirrorgeomlist(glist,'yz')
    glist = glist + mirrorgeomlist(glist,'xz')
    return glist

if __name__ == "__main__":

    from ezdxf_drawable import *
    #set up DXF rendering
    dd=ezdxfDraw()
    filename="example8-out"
    print("Output file name is {}.dxf".format(filename))
    dd.saveas(filename)

    glist = mirrorArray()

    drawLegend(dd)
    dd.draw(glist)
    dd.display()
