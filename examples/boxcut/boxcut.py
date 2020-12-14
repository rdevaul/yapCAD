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

import math
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
box_length = mindim
box_width = mindim
box_height = mindim

# base filename for DXF generation
outfile = "boxout"

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

def topPoly(toplen=box_length,topwid=box_width):
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

def frontPoly(lrlen=box_length,lrht=box_height):
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

def leftPoly(lrwid=box_width,lrht=box_height):
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

    top = topPoly(box_length,box_width)
    front = frontPoly(box_length,box_height)
    left = leftPoly(box_width,box_height)
    
    polylist.append(top)
    polylist.append(front)
    polylist.append(left)

    return polylist

## Draw some documentary text using the specified drawable.
## Text in the OpenGL rendering isn't working very well yet.

def legend(d):

    x = (box_length/2)*1.2
    y = (box_width/2+box_height)
    
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
    d.draw_text("length: {}mm".format(box_length),
                point(x,y),
                attr={'height': 6.0})
    y -= 8
    d.draw_text("width: {}mm".format(box_width),
                point(x,y),
                attr={'height': 6.0})
    y -= 8
    d.draw_text("height: {}mm".format(box_height),
                point(x,y),
                attr={'height': 6.0})


def twoDDraw(polys,dd,zoff=0.0):

    dd.layer = 'DOCUMENTATION'
    dd.linecolor = 'yellow'
    legend(dd)

    # calculate some spacing intervals
    lenhei = (box_length+box_height)/2.0
    widhei = (box_width+box_height)/2.0
    
    for i in range(len(polys)):
        p = polys[i]
        d1 = []
        d2 = []
        d1s = ""
        d2s = ""
        if i == 0:
            d1 = point(0,0,zoff) # location of bottom
            d2 = point(lenhei*1.1*2,0,zoff) # location of top
            d1s = "bottom"
            d2s = "top"
        elif i == 1: 
            d1 = point(0,-widhei*1.1,zoff) # location of front
            d2 = point(0,widhei*1.1,zoff) # location of back
            d1s = "front"
            d2s = "back"
        else: # i == 2
            d1 = point(-lenhei*1.1,0,zoff) #location of left
            d2 = point(lenhei*1.1,0,zoff) #location of right
            d1s = "left"
            d2s = "right"
            
        dd.layer = 'PATHS'
        dd.linecolor = 'white'

        g = p.geom()
        g1 = translate(g,d1)
        g2 = translate(g,d2)
        
        dd.draw(g1)
        dd.draw(g2)

        dd.layer = 'DOCUMENTATION'
        dd.linecolor = 'yellow'
        
        # att = {'style': 'LiberationMono','height': 10.0}
        att = {'height': 6.0, 'anchor_x': 'center'}
        
        dd.draw_text(d1s,d1, align='CENTER', attr=att)
        dd.draw_text(d2s,d2, align='CENTER', attr=att)
    
    
    if renderDebug:
        dd.linecolor='aqua'
        dd.draw(debugGeomList)

## So, yapCAD currently doesn't support rotating arcs out of the XY
## plane. That sucks for our visualization.  To hack around this, I'm
## creating a couple of utility functions that will replace arcs in
## geometry lists with sampled polyline representations.

# utility function to convert an arc to a polyline
def arc2poly(c, res=10.0):
    # resolution of arc sampling in degrees
    points = []
    p = c[0]
    r = c[1][0]
    start = c[1][1]
    end = c[1][2]
    if not (start==0 and end==360):
        start = start%360.0
        end = end % 360.0
        if end < start:
            end = end + 360
        theta = start*pi2/360.0
        points.append(add(p,[math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]))
        for a in range(round(start),round(end),round(res)):
            theta = a*pi2/360.0
            pp = [math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]
            pp = add(pp,p)
            points.append(pp)
        theta = end*pi2/360.0
        points.append(add(p,[math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]))
    return points

# utility function to replace arcs in a geometry list with sampled
# polyline represntations
def glistArc2poly(glist):
    ngl = []
    for e in glist:
        if isarc(e):
            ngl.append(arc2poly(e))
        else:
            ngl.append(e)
    return ngl

# take an XY poly, extrude in z, return geometry list. NOTE: this does
# not produce a surface representation, only a wireframe of lines as a
# geometry list

def extrudeZ(poly,thick):
    geom = poly.geom()
    geom = glistArc2poly(geom)
    # bottom = reverseGeomList(geom)
    bottom = geom
    bottom = translate(bottom,point(0,0,-thick/2))
    top = geom
    top = translate(top,point(0,0,thick/2))

    def glist2struts(glist):
        struts=[]
        upr = point(0,0,thick)
        for e in glist:
            strut=[]
            if ispoint(e):
                strut = line(e,add(e,upr))
                struts.append(strut)
            elif isline(e):
                strut = line(e[0],add(e[0],upr))
                struts.append(strut)
            elif isarc(e):
                p = samplearc(e,0.0)
                strut = line(p,add(p,upr))
                struts.append(strut)
            elif ispoly(e) or isgeomlist(g):
                struts = struts + glist2struts(e)
            else:
                raise ValueError('bad thing passed to glist2struts')
        return struts
    
    strts = glist2struts(bottom)
    return bottom + top + strts
    
    
def makeOglGeom():
    ## regenerate geometry with zero kerf correction
    global kerf
    kerf = 0.0
    polys = makePolys()

    parts = []
    for poly in polys:
        part = extrudeZ(poly,thick)
        parts.append(part)

    model = []

    for i in range(len(parts)):
        part = parts[i]
        comp1 = []
        comp2 = []
        if i == 0: # top and bottom
            comp1 = translate(part,point(0,0,-(box_height-thick/2)/2))
            comp2 = translate(part,point(0,0,(box_height-thick/2)/2))
        elif i == 1: # front and back
            part = rotate(part,90.0,axis=point(1,0,0))
            comp1 = translate(part,point(0,-(box_width-thick/2)/2,0))
            comp2 = translate(part,point(0,(box_width-thick/2)/2,0))
        else: # i = 2, left and right
            part = rotate(part,90.0,axis=point(0,1,0))
            comp1 = translate(part,point(-(box_length-thick/2)/2,0,0))
            comp2 = translate(part,point((box_length-thick/2)/2,0,0))
        model.append(comp1)
        model.append(comp2)

    return model
        
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

    box_width = args.width*dimmul
    box_length = args.length*dimmul
    box_height = args.height*dimmul

    if box_width < mindim or box_length < mindim or box_height < mindim:
        raise ValueError('minimum length of {}mm required for each dimension'.format(mindim))

    outname = args.outname
    if outname == "":
        raise ValueError('non-zero length file name base required')
    
    if args.debug:
        renderDebug=True

    if args.thick > 0:
        thick = args.thick

    if args.kerf >= 0:
        kerf = args.kerf

    if args.gl:
        ogl = True

    polys = makePolys()
    dd = ezdxfDraw()
    dd.filename = outfile

    twoDDraw(polys,dd)
    dd.display()

    if ogl:
        from yapcad.pyglet_drawable import *
        dd2 = pygletDraw()
        # magnification factor for text
        dd2.magnify=1.5
        # compute the camera distance assuming a 60 degree (pi/6) FOV
        maxd = max((box_length*2+box_height*2)*1.1,
                   (box_width+ box_height*2)*1.1)
        dist = ((maxd/2) / math.sin(pi/6))*math.cos(pi/6)
        dd2.cameradist = dist-box_height*2

        
        twoDDraw(polys,dd2,-box_height*2)

        oglGeom = makeOglGeom()
        #print("oglGeom: ",oglGeom)

        #oglGeom = translate(oglGeom,point(0,0,box_height))

        dd2.layer='PATHS'
        for i in range(len(oglGeom)):
            if i == 0 or i == 1:
                dd2.linecolor = 'silver'
            elif i == 2 or i == 3:
                dd2.linecolor = 'aqua'
            else:
                dd2.linecolor = 'fuchsia'
                
            dd2.draw(oglGeom[i])
        dd2.display()
    print("done")

