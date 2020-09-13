## yapCAD 3D geometry example

print("example10.py -- yapCAD inside-outside testing demonstration")
print("""
In this example, we create different types of yapCAD geometry and
randomly gemerate test ponts.  Test points that fall inside are
rendered as aqua crosses, points that fall outside are rendered as red
crosses.
""")

from yapcad.geom import *
import random

def randomPoints(bbox,numpoints):
    points = []
    minx = bbox[0][0]
    maxx = bbox[1][0]
    miny = bbox[0][1]
    maxy = bbox[1][1]
    rangex = maxx-minx
    rangey = maxy-miny
    for i in range(numpoints):
        points.append(point(random.random()*rangex+minx,
                            random.random()*rangey+miny))
    return points

def pointInBox(bbox,r):
    minx = bbox[0][0]+r
    maxx = bbox[1][0]-r
    miny = bbox[0][1]+r
    maxy = bbox[1][1]-r
    rangex = maxx-minx
    rangey = maxy-miny
    x = point(random.random()*rangex+minx,
              random.random()*rangey+miny)
    return x


def randomArc(bbox,minr=0.0,maxr=10.0,circle=False):
    radr = maxr-minr
    r = random.random()*radr+minr
    start = 0
    end = 360
    if not circle:
        start = random.random()*360.0
        end = start + random.random()*360.0

    x = pointInBox(bbox,r)
    return arc(x,r,start,end)

def randomPoly(bbox,numpoints=10,minr = 1.0,maxr = 10.0):
    angles = []
    rads = []
    ang = 0.0
    rr = maxr-minr
    for i in range(numpoints):
        a = random.random()
        angles.append(ang)
        ang = ang + a
        # print("a: ",a," ang: ",ang)
        rads.append(random.random()*rr+minr)

    sf = pi2/ang

    points = []
    x = pointInBox(bbox,maxr)
    
    for i in range(numpoints):
        p = [cos(angles[i]*sf)*rads[i],
             sin(angles[i]*sf)*rads[i],0,1]
        points.append(add(p,x))

    return points + [ points[0] ]

def randomGeometry(bbox, numarcs, numcircles, numpolys):
    geom=[]
    for i in range(numarcs):
        geom.append(randomArc(bbox,minr = 5,maxr=10))

    for i in range(numcircles):
        geom.append(randomArc(bbox,minr = 3, maxr=10,circle=True))

    for i in range(numpolys):
        geom.append(randomPoly(bbox,numpoints=round(random.random()*7+3),
                               maxr =10))
    return geom

def drawGeom(dd,geom):
    for g in geom:
        if isarc(g) and not iscircle(g):
            ## draw the piza slice
            dd.draw(g)
            dd.draw(line(sample(g,1.0),
                         g[0]))
            dd.draw(line(g[0],
                         sample(g,0.0)))
        else:
            dd.draw(g)

def testPoints(points,geom):
    inpts=[]
    outpts=[]

    for p in points:
        ins = False
        for g in geom:
            if inside(g,p):
                inpts.append(p)
                ins = True
                continue
        if not ins:
            outpts.append(p)
    return inpts, outpts

if __name__ == "__main__":
    from yapcad.pyglet_drawable import *
    dd=pygletDraw()
    dd.set_linecolor('white')

    bbox = line([-60,-60,0,1],[60,60,0,1])

    tps = randomPoints(bbox,2000)

    #glist = randomGeometry(bbox,20,20,10)
    glist = randomGeometry(bbox,20,20,20)

    inpts,outpts = testPoints(tps,glist)

    drawGeom(dd,glist)

    dd.polystyle='points'

    dd.set_linecolor('aqua')
    dd.pointstyle = 'x'
    dd.draw(inpts)

    dd.set_linecolor('red')
    dd.pointstyle = 'x'
    dd.draw(outpts)
    dd.display()
