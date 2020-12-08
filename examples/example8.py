## yapCAD poly() intersecton and drawing examples

from yapcad.geom import *
from yapcad.poly import *
import examples.example10 as example10
import random

def drawLegend(d):
    ## Put some documentary text on the drawing
    d.layer = 'DOCUMENTATION'

    att = {'style': 'OpenSans-Bold',
           'height': 2.5}
    
    d.draw_text("yapCAD", point(5,15), attr = att)

    d.draw_text("example8.py", point(5,11), attr=att)
    d.draw_text("Polygon() flowers, mirrored geometry",
                point(5,7),attr=att)
    d.layer = False # back to default layer


def flower(center = point(0,0),
           petals = 10,
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
    flwr = flower(pnt,returnPoly=True)
    glist = flwr.geom()
    bb = bbox(glist)
    flwr2 = deepcopy(flwr)
    flwr2.grow(1.0)
    glist.append(flwr2.geom())

    ranp = example10.randomPoints(bb,500)

    for p in ranp:
        if flwr.isinside(p):
            glist.append(arc(p,0.4))
                         
    ply = [point(bb[0]), point(bb[1][0],bb[0][1]),
           point(bb[1]), point(bb[0][0],bb[1][1])]
    glist = glist + ply
    glist = glist + mirror(glist,'yz')
    glist = glist + mirror(glist,'xz')
    return glist

if __name__ == "__main__":

    from yapcad.ezdxf_drawable import *
    #set up DXF rendering

    print("example8.py -- yapCAD computational geometry and DXF drawing example")
    print("""
This is a demo of the Polygon() and mirror capabilities of yapCAD,
making and then mirroring random "flowers."
    """)
    dd=ezdxfDraw()
    filename="example8-out"
    print("Output file name is {}.dxf".format(filename))
    dd.filename = filename

    glist = mirrorArray()

    def mydrawer(gl):
        g1=[]
        g2=[]
        for g in gl:
            if iscircle(g):
                g1.append(g)
            else:
                g2.append(g)
        dd.linecolor = 'white'
        dd.draw(g2)
        dd.linecolor = 'aqua'
        dd.draw(g1)

    drawLegend(dd)
    mydrawer(glist)
    dd.display()
