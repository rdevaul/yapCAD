## multi-rendering-back-end drawing example for yapCAD
print("example1.py -- yapCAD DXF and OpenGL drawing example")

from yapcad.ezdxf_drawable import *
from yapcad.pyglet_drawable import *
from yapcad.geom import *

#set up openGL rendering
def setupGL():
    dGl =pygletDraw()
    dGl.magnify = 1.0
    dGl.linecolor = 'white'
    dGl.cameradist = 25
    return dGl

#set up DXF rendering
def setupDXF():
    d=ezdxfDraw()
    filename="example1-out"
    print("\nOutput file name is {}.dxf".format(filename))
    d.filename = filename
    return d

## Draw some documentary text using the specified drawable.
## Text in the OpenGL rendering isn't working very well yet.

def legend(d):

    d.draw_text("yapCAD", point(5,15),\
                attr={'style': 'OpenSans-Bold', # style for ezdxf
                      'font_name': 'OpenSans', # style for pyglet
                      'bold': True, # style for pyglet
                      'height': 2.0})
    d.draw_text("example1.py",
                point(5,12))

    d.draw_text("drawing primitives",
                point(5,10))

    
## make geometry for drawing, and return it as a geometry list
def geometry():
    p = point(10,8)
    l = line(point(-5,15),
             point(10,-5))
    a = arc(arc(point(1,3),6,45,135))

    return [p,l,a]


## Do the rendering with the specified drawable
## Note: if we just wanted to render all of this the same way,
## we could just call d.draw(glist)

def drawGlist(glist,d):

    ## extract the point, line, and arc from the list
    p = glist[0]
    l = glist[1]
    a = glist[2]

    ## in this example we show three different ways of specifying the
    ## line color used for rendering.  By color index (based on the
    ## AutoCAD standard color index table), by standard HTML color
    ## name, and by specifying the RGB tripple.  All of these methods
    ## should work regardless of the drawable back-end used.
    
    ## render the point
    d.pointstyle = 'xo' # set the point rendering style
    d.linecolor = 1 # set color to red (DXF index color)
    d.draw(p)

    ## render the line in white
    d.linecolor = 'white' #set color to white, by name
    d.draw(l)

    ## render the arc in aqua, by specifying the RGB tripple.
    ## NOTE:
    ## only RGB tripples that correspond to AutoCAD color table
    ## entries will produce the desired results in DXF rendering at
    ## present.
    
    d.linecolor = [0,255,255] # corresponds to 'aqua'
    d.draw(a)

if __name__ == "__main__":

    ## setup DXF rendering
    d= setupDXF()

    ## setup OpenGL rendering
    dGl = setupGL()

    ## make the geometry
    geomlist = geometry()

    ## draw the geometry in DXF
    drawGlist(geomlist,d)

    ## add a dawing legend on the DXF drawing, in the DOCUMENTATION
    ## layer
    
    d.layer = 'DOCUMENTATION'
    d.linecolor = 256 ## set layer default color
    legend(d)

    ## add a drawing legend in the OpenGL rendering, setting the
    ## color eplicitly to yellow
    
    dGl.linecolor = 'yellow'
    legend(dGl)
    
    ## draw the geometry in OpenGL
    drawGlist(geomlist,dGl)

    ## write out the DXF file as example1-out.dxf
    d.display()

    ## create the interactive OpenGL rendering -- do this last
    dGl.display()
    
