## yapCAD boolean operations geometry example

from yapcad.geom import *
from yapcad.octtree import *
import examples.example10 as example10
from yapcad.pyglet_drawable import *

docGeomList=[]

def draw(dd):
    foo2 = NTree()
    # make 1000 lines in space, keep track of which fall inside a
    # sub-region, then try to retrieve those with our tree,

    bigbox = [[-15,-15,-15,1],[15,15,15,1]]
    deltabox = [[-10,-10,-10,1],[10,10,10,1]]
    subbox = [[10,-10,5,1],[20,0,15,1]]
    glist=[]
    sublist=[]
        
    dim = 400
    plist1 = example10.randomPoints(bigbox,dim)
    plist2 = example10.randomPoints(deltabox,dim)
        
    for i in range(dim):
        l = line(plist1[i],add(plist1[i],plist2[i]))
        lbx = bbox(l)
        if boxoverlap2(lbx,subbox):
            sublist.append(l)
        else:
            glist.append(l)
        foo2.addElement(l)

    dd.linecolor = [127,79,63]
    dd.draw(glist)
    dd.linecolor = 'white'
    dd.draw_bbox(subbox,dim3=True)
    dd.linecolor = 'aqua'
    dd.draw(sublist)
    dd.linecolor = 'red'
    for l in sublist:
        dd.draw_bbox(l,dim3=True)

    foo2.updateTree()
    print ("foo2 depth: ",foo2.depth)
    slist2 = foo2.getElements(subbox)
    assert len(slist2) == len(sublist)
    for l in slist2:
        assert l in sublist
    

def drawLegend(d):
    ## Put some documentary text on the drawing
    d.layer = 'DOCUMENTATION'

    att = {'style': 'OpenSans-Bold',
           'height': 2.5}
    
    d.draw_text("yapCAD", point(5,15), attr = att)

    d.draw_text("example13.py", point(5,11), attr=att)
    d.draw_text("octree experiment",
                point(5,7),attr=att)
    d.layer = False # back to default layer


if __name__ == "__main__":
    import sys
    renderOgl = False
    filename="example13-out"
    print("example13.py -- yapCAD boolean operations demonstration")

    oglarg= ("pyglet","opengl","OpenGL")
    dxfarg= ("ezdxf","dxf","DXF")
    if len(sys.argv) > 1:
        if sys.argv[1] in oglarg:
            renderOgl=True
        elif sys.argv[1] in dxfarg:
            renderOgl=False
        else:
            print("This is an experiment with octrees")
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
    
    draw(dd)
    dd.display()
    print("done")
