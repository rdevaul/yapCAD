## yapCAD boolean operations geometry example

from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *
import examples.example8 as example8

docGeomList=[]

def makeJigsawPiece(width=100,height=50,diam=30,frac=0.80,upper=False,poly=False):
    rect = makeRoundRect(width,height,5.0)
    rect = rect.translate(point(0,-height/2),poly=True)
    circ1 = makeCircle(point(-width/4,(diam/2)*frac),diam/2)
    #circ2 = makeCircle(point(width/4,-(diam/2)*frac),diam/2)
    circ2 = makeRoundRect(diam,diam,20.0,center=point(width/4,-(diam/2)*frac))

    docGeomList.append(circ1.geom())
    docGeomList.append(circ2.geom())
    bln = Boolean('union',[circ1,rect])
    bln = Boolean('difference',[bln,circ2])
    docGeomList.append(bln.geom())
    if poly:
        return bln
    return bln.geom()

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

    logopoly = example8.flower(returnPoly=True)
    jigsaw=makeJigsawPiece(poly=True)
    bigc= makeCircle(point(0,0,0),20)
    nb = Boolean('union',[logopoly,bigc])
    nb = nb.translate(point(0,-50),poly=True)
    docGeomList.append(nb.geom())

    nb = Boolean('difference',[jigsaw,nb])
    
    gl = nb.geom()
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

    if not renderOgl:
        dd.layer = 'DOCUMENTATION'
        dd.linecolor = False
    else:
        dd.linecolor = 'yellow'
    dd.draw(docGeomList)

    if not renderOgl:
        dd.layer = 'PATHS'
        
    dd.linecolor = 'white'
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
            print("Boolean() instances. ")
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
        dd.cameradist=170.0
        
    else:
        print("DXF rendering selected")
        from yapcad.ezdxf_drawable import *
        #set up DXF rendering
        dd=ezdxfDraw()
        dd.filename = filename
        drawLegend(dd)
    print("rendering...")
    
    testAndDraw(dd)
    print("done")
