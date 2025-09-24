## yapCAD boolean operations and surface creation geometry example

from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *
from yapcad.geom_util import *
from yapcad.geom3d import *

import examples.example8 as example8

docGeomList=[]

def makeJigsawPiece(width=100,height=50,diam=30,frac=0.80,upper=False,poly=False):
    rect = RoundRect(width,height,5.0)
    rect.translate(point(0,-height/2))
    circ1 = Circle(point(-width/4,(diam/2)*frac),diam/2)
    #circ2 = Circle(point(width/4,-(diam/2)*frac),diam/2)
    circ2 = RoundRect(diam,diam,20.0,center=point(width/4,0))

    docGeomList.append(circ1.geom)
    docGeomList.append(circ2.geom)
    bln = Boolean('union',[circ1,rect])
    bln = Boolean('difference',[bln,circ2])
    docGeomList.append(bln.geom)
    if poly:
        return bln
    return bln.geom

def drawLegend(d):
    ## Put some documentary text on the drawing
    d.layer = 'DOCUMENTATION'

    att = {'style': 'OpenSans-Bold',
           'height': 2.5}
    
    d.draw_text("yapCAD", point(5,15), attr = att)

    d.draw_text("example12.py", point(5,11), attr=att)
    d.draw_text("complex boolean operations on Polygon() instances",
                point(5,7),attr=att)
    d.layer = False # back to default layer


def geometry():

    logopoly = example8.flower(maxRadius=30,
                               returnPoly=True)
    jigsaw=makeJigsawPiece(poly=True)
    bigc= Circle(point(0,0,0),20)
    nb = Boolean('union',[logopoly,bigc])
    nb.translate(point(0,-50))
    docGeomList.append(nb.geom)

    nb = Boolean('difference',[jigsaw,nb])
    
    gl = nb.geom
    gl2 = []

    for i in range(100):
        u = i/99.0
        rad = 1.0 + u*2
        
        p = nb.sample(u)
        c = arc(p,rad)
        gl2.append(c)
        

    return gl,gl2


def testAndDraw(dd):
    geom,geom2=geometry()

    if not renderOgl:
        dd.layer = 'DRILLS'
        dd.linecolor = False

    else:
        dd.linecolor = 'aqua'
    dd.draw(geom2)
    if renderOgl:
        dd.linecolor = 'red'
        dd.draw(translate(geom,point(0,0,5)))
    else:
        dd.draw(geom)

    if not renderOgl:
        dd.layer = 'DOCUMENTATION'
        dd.linecolor = False
    else:
        dd.linecolor = 'red'
    dd.draw(docGeomList)

    if not renderOgl:
        dd.layer = 'PATHS'
        
    dd.linecolor = 'white'
    if renderOgl:
        ply = geomlist2poly(translate(geom,point(0,0,-0.1)))
        surf,bnd = poly2surfaceXY(ply)
        dd.draw_surface(surf)
    
    geom = translate(geom,point(0,0,0.1))
    dd.draw(geom)

    dd.display()
    

if __name__ == "__main__":
    import sys
    renderOgl = False
    filename="example12-out"
    print("example12.py -- yapCAD boolean operations demonstration")

    oglarg= ("pyglet","opengl","OpenGL")
    dxfarg= ("ezdxf","dxf","DXF")
    if len(sys.argv) > 1:
        if sys.argv[1] in oglarg:
            renderOgl=True
        elif sys.argv[1] in dxfarg:
            renderOgl=False
        else:
            print("This is an example of a complex boolean operation on Polygon() and")
            print("Boolean() instances. When rendering with OpenGL, we also create a")
            print("two-dimensional surface from the resulting figure.")
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
            dd.cameradist=170.0
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
        drawLegend(dd)
        renderOgl = False
    print("rendering...")
    
    testAndDraw(dd)
    print("done")
