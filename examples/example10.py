## yapCAD 3D inside-outside testing example
## Copyright (c) 2020 Richard DeVaul

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
In this example, we create different types of yapCAD geometry and
randomly gemerate test points.  Test points that fall inside are
rendered as aqua crosses, points that fall outside are rendered as red
crosses.
"""

from yapcad.geom import *
from yapcad.geom_util import *
import yapcad.poly as poly
import random


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

def testPoints(points,geom,testElements=True):
    inpts=[]
    outpts=[]

    for p in points:
        ins = False
        if testElements:
            for g in geom:
                if isinsideXY(g,p):
                    inpts.append(p)
                    ins = True
                    continue
        else:
            if isinstance(geom,poly.Polygon):
                if geom.isinsideXY(p):
                    inpts.append(p)
                    ins = True
            else:
                if isinsideXY(geom,p):
                    inpts.append(p)
                    ins = True
        if not ins:
            outpts.append(p)
    return inpts, outpts

def testAndDraw(dd):
    dd.linecolor = 'white'

    ## this is the bounding box for our test area
    bbox = line([-60,-60,0,1],[60,60,0,1])

    ## make 2,000 points
    tps = randomPoints(bbox,2000)

    ## make 20 random arcs, 20 random circles, and 20 random polys.
    glist = randomGeometry(bbox,20,20,20)

    ## test the list of points, sort in to inside points and outside
    ## points
    inpts,outpts = testPoints(tps,glist)

    ## draw the test geometry -- we are using a special draw function
    ## that will render arcs as "pizza slices."
    drawGeom(dd,glist)

    ## Draw the points.
    dd.polystyle='points'

    dd.linecolor = 'aqua'
    dd.pointstyle = 'x'
    dd.draw(inpts)

    dd.linecolor = 'red'
    dd.pointstyle = 'o'
    dd.draw(outpts)
    dd.display()
    
if __name__ == "__main__":
    import sys
    print("example10.py -- yapCAD inside-outside testing demonstration")

    renderOgl = False
    filename="example10-out"
    oglarg= ("pyglet","opengl","OpenGL")
    dxfarg= ("ezdxf","dxf","DXF")
    if len(sys.argv) > 1:
        if sys.argv[1] in oglarg:
            renderOgl=True
        elif sys.argv[1] in dxfarg:
            renderOgl=False
        else:
            print("""
 In this example, we create different types of yapCAD geometry and
 randomly gemerate test points.  Test points that fall inside are
 rendered as aqua crosses, points that fall outside are rendered as red
 crosses.
            """)

            print("syntax: $ python3 {} <rendertype> [filename]".format(sys.argv[0]))
            print("    where <rendertype> is one of {} for OpenGL".format(oglarg))
            print("    or one of {} for DXF".format(dxfarg))
            print("    For DXF, you can optionally specify [filename].dxf as the output file")
            quit()
    if len(sys.argv) > 2 and renderOgl==False:
        filename = sys.argv[2]+".dxf"
    requestedOgl = renderOgl
    dd = None
    if renderOgl:
        print("OpenGL rendering selected")
        try:
            from yapcad.pyglet_drawable import *
            dd=pygletDraw()
        except RuntimeError as err:
            print("\nSkipping OpenGL rendering: {}".format(err))
            renderOgl = False
        except Exception as err:
            print("\nSkipping OpenGL rendering due to unexpected error: {}".format(err))
            renderOgl = False

    if not renderOgl or dd is None:
        if requestedOgl and dd is None:
            print("Falling back to DXF rendering")
        else:
            print("DXF rendering selected")
        from yapcad.ezdxf_drawable import *
        #set up DXF rendering
        dd=ezdxfDraw()
        dd.filename = filename
        renderOgl = False
    print("rendering...")
    testAndDraw(dd)
    print("done")
