## boxcut procedural design example for yapCAD
## Born on 12 December, 2020
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

"""procedural design example that generates squeeze-fit box designs from
input parameters and produces DXF files as output.

This is intended as a demonstration design system -- more work would
be needed to generalize this to designs bigger than 200mm or so on a
side.

"""

from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *
from yapcad.ezdxf_drawable import *

## global configuration variables

# minimum box dimension, mm
# NOTE: this is dictated by the tab structure
mindim = 65.0

#tab length, mm
tlen = 6.0

## global state variables, set in __main__

# box dimensions, nominal default values.  Internally, all dimensions
# are in mm.
length = 10.0
width = 10.0
height = 10.0

# base filename for DXF generation
outfile = "boxout"

# flag for generating separate drawings rather than integrated
separate = False

# material thickness, mm
thick = 3.175

# kerf compensation, mm
kerf = 0.397

# do we do openGL?
ogl = False

# a global geometry list we use for rendering debug-useful geometry
debugGeomList=[]

renderDebug=False

# define the bottom and top pieces by outline.  These share the same
# edge-join geometry:

#      +----__--__------__--__----+ 
#      |                          | 
#      |                      ^   | 
#       |  Bottom and Top     |  |  
#      |                      |   | 
#       |  Tab pattern is     |  |    
#      |   the same on      width | 
#      |   all edges.         |   |
#       |                     |  |  
#      | <------ length ------+-> |
#       |                     |  |
#      |                      v   |
#      |    __  __      __  __    |
#      +----  --  ------  --  ----+

def topPoly(toplen=length,topwid=width):
    # start with a rounded rectangle the size of the face
    face = makeRoundRect(toplen+kerf*2,topwid+kerf*2,2.0)

    debugGeomList.append(face.geom())
    # make trimming polygons
    trimtab = makeRect(tlen-kerf*2,(thick-kerf)*2)
    trimtab2 = makeRect((thick-kerf)*2,tlen-kerf*2)

    tl2 = toplen/2.0
    tw2 = topwid/2.0
    t1off = tlen*2.5
    t2off = tlen*4.5

    # duplicate and position to trim the upper, right corner
    t1a = trimtab.translate(point(tl2-t1off,
                                  tw2),poly=True)
    t1b = trimtab.translate(point(tl2-t2off,
                                  tw2),poly=True)
    t2a = trimtab2.translate(point(tl2,
                                   tw2-t1off),poly=True)
    t2b = trimtab2.translate(point(tl2,
                                   tw2-t2off),poly=True)

    debugGeomList.append(t1a.geom())
    debugGeomList.append(t1b.geom())
    debugGeomList.append(t2a.geom())
    debugGeomList.append(t2b.geom())
    
    # perform boolean operations for each corner -- use mirror to
    # rotate the trimming bodies around the rectangular form.

    for i in range(4):
        if i == 0:
            pass # we are in correct position for first trim
        elif i == 1 or i == 3:
            t1a = t1a.mirror("yz",poly=True)
            t1b = t1b.mirror("yz",poly=True)
            t2a = t2a.mirror("yz",poly=True)
            t2b = t2b.mirror("yz",poly=True)
        else: # i == 2
            t1a = t1a.mirror("xz",poly=True)
            t1b = t1b.mirror("xz",poly=True)
            t2a = t2a.mirror("xz",poly=True)
            t2b = t2b.mirror("xz",poly=True)
            
        face = Boolean('difference',[face,t1a])
        face = Boolean('difference',[face,t1b])
        face = Boolean('difference',[face,t2a])
        face = Boolean('difference',[face,t2b])

    return face

# define the front and back pieces by outline.  These share the same
# edge-join geometry:

#           __  __      __  __      
#      +----  --  ------  -- ^----+ 
#      |                     |    | 
#       |  Front and Back    |   |  
#      |                     |    | 
#       |  Tab pattern is    |   |    
#      |   the same on     height | 
#      |   opposite edges    |    |
#       |                    |   |  
#      | <------ length -----+--> |
#       |                    |   |
#      |                     v    |
#      +----__--__------__--__----+
#                                  

def frontPoly(lrlen=length,lrht=height):
    # start with a rounded rectangle the size of the face, less the
    # projecting tabs
    face = makeRoundRect(lrlen+kerf*2,lrht+(kerf-thick)*2,2.0)

    debugGeomList.append(face.geom())
    # make tab and trim polygons
    tab = makeRect(tlen+kerf*2,(thick+kerf)*2)
    trimtab = makeRect((thick+kerf)*2,tlen+kerf*2)

    tl2 = lrlen/2.0
    tw2 = lrht/2.0
    t1off = tlen*2.5
    t2off = tlen*4.5

    # duplicate and position to trim the upper, right corner
    t1a = tab.translate(point(tl2-t1off,
                              tw2-thick),poly=True)
    t1b = tab.translate(point(tl2-t2off,
                              tw2-thick),poly=True)
    t2a = trimtab.translate(point(tl2,
                                  tw2-t1off),poly=True)
    t2b = trimtab.translate(point(tl2,
                                  tw2-t2off),poly=True)

    debugGeomList.append(t1a.geom())
    debugGeomList.append(t1b.geom())
    debugGeomList.append(t2a.geom())
    debugGeomList.append(t2b.geom())
    
    # perform boolean operations for each corner -- use mirror to
    # rotate the trimming bodies around the rectangular form.

    for i in range(4):
        if i == 0:
            pass # we are in correct position for first trim
        elif i == 1 or i == 3:
            t1a = t1a.mirror("yz",poly=True)
            t1b = t1b.mirror("yz",poly=True)
            t2a = t2a.mirror("yz",poly=True)
            t2b = t2b.mirror("yz",poly=True)
        else: # i == 2
            t1a = t1a.mirror("xz",poly=True)
            t1b = t1b.mirror("xz",poly=True)
            t2a = t2a.mirror("xz",poly=True)
            t2b = t2b.mirror("xz",poly=True)
            
        face = Boolean('union',[face,t1a])
        face = Boolean('union',[face,t1b])
        face = Boolean('difference',[face,t2a])
        face = Boolean('difference',[face,t2b])

    return face

# define the left and right pieces by outline.  These share the same
# edge-join geometry:

#           __  __      __  __      
#       +---  --  ------  -- ^---+  
#       |                    |   |  
#      |   Left and Right    |    | 
#       |                    |   |  
#      |   Tab pattern is    |    |   
#       |  the same on     height|  
#       |  opposite edges    |   | 
#      |                     |    | 
#       |                    |   | 
#      |<-------- width -----+--->|
#       |                    v   | 
#       +---__--__------__--__---+ 
#                                  

def leftPoly(lrwid=width,lrht=height):
    # start with a rounded rectangle the size of the face, less the
    # projecting tabs
    face = makeRoundRect(lrwid+(kerf-thick)*2,lrht+(kerf-thick)*2,2.0)

    debugGeomList.append(face.geom())
    # make tab and trim polygons
    tab1 = makeRect(tlen+kerf*2,(thick+kerf)*2)
    tab2 = makeRect((thick+kerf)*2,tlen+kerf*2)

    tl2 = lrwid/2.0
    tw2 = lrht/2.0
    t1off = tlen*2.5
    t2off = tlen*4.5

    # duplicate and position to add tabs to the upper, right corner
    t1a = tab1.translate(point(tl2-t1off,
                              tw2-thick),poly=True)
    t1b = tab1.translate(point(tl2-t2off,
                              tw2-thick),poly=True)
    t2a = tab2.translate(point(tl2-thick,
                               tw2-t1off),poly=True)
    t2b = tab2.translate(point(tl2-thick,
                               tw2-t2off),poly=True)

    # add these to our debugging geometry list for later
    # visualization, if desired
    debugGeomList.append(t1a.geom())
    debugGeomList.append(t1b.geom())
    debugGeomList.append(t2a.geom())
    debugGeomList.append(t2b.geom())
    
    # perform boolean operations for each corner -- use mirror to
    # rotate the trimming bodies around the rectangular form.

    for i in range(4):
        if i == 0:
            pass # we are in correct position for first corner
        elif i == 1 or i == 3:
            t1a = t1a.mirror("yz",poly=True)
            t1b = t1b.mirror("yz",poly=True)
            t2a = t2a.mirror("yz",poly=True)
            t2b = t2b.mirror("yz",poly=True)
        else: # i == 2
            t1a = t1a.mirror("xz",poly=True)
            t1b = t1b.mirror("xz",poly=True)
            t2a = t2a.mirror("xz",poly=True)
            t2b = t2b.mirror("xz",poly=True)
            
        face = Boolean('union',[face,t1a])
        face = Boolean('union',[face,t1b])
        face = Boolean('union',[face,t2a])
        face = Boolean('union',[face,t2b])

    return face.rotate(90,poly=True)

def makePolys():
    polylist = []

    top = topPoly(length,width)
    front = frontPoly(length,height)
    left = leftPoly(width,height)
    
    polylist.append(top)
    polylist.append(front)
    polylist.append(left)

    return polylist

## Draw some documentary text using the specified drawable.
## Text in the OpenGL rendering isn't working very well yet.

def legend(d):

    x = (length/2)*1.2
    y = (width/2+height)
    
    d.draw_text("yapCAD", point(x,y),\
                attr={'style': 'OpenSans-Bold', # style for ezdxf
                      'font_name': 'OpenSans', # style for pyglet
                      'bold': True, # style for pyglet
                      'height': 20.0})
    y -= 25
    d.draw_text("boxcut.py",
                point(x,y),
                attr={'height': 15.0})

    y -= 20
    d.draw_text("length: {}mm".format(length),
                point(x,y),
                attr={'height': 6.0})
    y -= 8
    d.draw_text("width: {}mm".format(width),
                point(x,y),
                attr={'height': 6.0})
    y -= 8
    d.draw_text("height: {}mm".format(height),
                point(x,y),
                attr={'height': 6.0})


def dxfDraw(polys):
    dd = ezdxfDraw()
    dd.filename = outfile

    dd.layer = 'DOCUMENTATION'
    legend(dd)
    if renderDebug:
        dd.linecolor='aqua'
        dd.draw(debugGeomList)

    dd.layer = 'PATHS'
    dd.linecolor = 'white'

    # calculate some spacing intervals
    lenhei = (length+height)/2.0
    widhei = (width+height)/2.0
    
    for i in range(len(polys)):
        p = polys[i]
        d1 = []
        d2 = []
        if i == 0:
            d1 = point(0,0) # location of bottom 
            d2 = point(lenhei*1.1*2,0) # location of top
        elif i == 1: 
            d1 = point(0,-widhei*1.1) # location of front
            d2 = point(0,widhei*1.1) # location of back
        else: # i == 2
            d1 = point(-lenhei*1.1,0) #location of left
            d2 = point(lenhei*1.1,0) #location of right
            
        dd.draw(p.translate(d1))
        dd.draw(p.translate(d2))

    dd.display()

def renderOgl(polys):
    print("draw openGL stuff here")
    

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="""
    procedural design example that generates squeeze-fit box designs from
    input parameters and produces DXF files as output.
    """)

    parser.add_argument("length",type=float, help="length of the box")
    parser.add_argument("width",type=float, help="width of the box")
    parser.add_argument("height",type=float, help="height of the box")
    parser.add_argument("-u","--units",type=str,choices=["mm","inch"],
                        default="mm",help="units for dimensions")
    parser.add_argument("-o","--outname",type=str,default="boxout",
                        help="file name base for output, default is boxout")
    parser.add_argument("-s","--separate",action="store_true",
                        help="make one DXF drawing per face rather than a combined tiled layout")
    parser.add_argument("-t","--thick",type=float,default=-1,
                        help="material thickness, default is 3.175 mm (1/4\")")
    parser.add_argument("-k","--kerf",type=float,default=-1,
                        help="kerf compensation, default is 0.397mm (1/64\")")
    parser.add_argument("-g","--gl",action="store_true",
                        help="do an interactive 3D visualization of the box")
    parser.add_argument("-d","--debug",action="store_true",
                        help="turn on rendering of debug geometry")
    args = parser.parse_args()
    print(args)

    dimmul = 1.0
    if args.units == "inch":
        dimmul = 25.4

    width = args.width*dimmul
    length = args.length*dimmul
    height = args.height*dimmul

    if width < mindim or length < mindim or height < mindim:
        raise ValueError('minimum length of {}mm required for each dimension'.format(mindim))

    outname = args.outname
    if outname == "":
        raise ValueError('non-zero length file name base required')
    
    if args.separate:
        separate=True

    if args.debug:
        renderDebug=True

    if args.thick > 0:
        thick = args.thick

    if args.kerf >= 0:
        kerf = args.kerf

    if args.gl:
        ogl = True

    polys = makePolys()
    dxfDraw(polys)

    if ogl:
        renderOgl(polys)
        
    print("done")

