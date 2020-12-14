## simple yapCAD framework for openGL drawing using pyglet
## package
## Copyright (c) 2020 Richard W. DeVaul
## Copyright (c) 2020 yapCAD contributors
## All rights reserved

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

import math
import pyglet
import pyglet.font as font
import pyglet.gl as gl
import pyglet.graphics as graphics

from yapcad.geom import *
import yapcad.drawable as drawable

## openGL utility functions
def vec(*args):
    return (gl.GLfloat * len(args))(*args)

## HTML document, instructions for user interaction
yapCAD_legend="""
<font face="OpenSans, Geneva, sans-serif" size="5" color="#f0f0ff"><b>yapCAD</b><br></font>
<font face="Verdana, Geneva, sans-serif" size="3" color="white"><b>3D drawing viewer</b></font><br>
<font face="Verdana, Geneva, sans-serif" size="1" color="white">
<br>
<b>up-arrow</b>: zoom in<br>
<b>down-arrow</b>: zoom out<br>
<b>l</b>: toggle lighting mode<br>
<b>p</b>: toggle ground plane<br>
<b>left-mouse drag</b>: rotate view<br>
<b>m</b>: toggle display of this message<br>
<b>ESC</b>: exit viewer<br>
</font>
"""

## class to provide openGL drawing functionality
class pygletDraw(drawable.Drawable):

    def window(self):
        window = []
        try:
            # Try and create a window with multisampling (antialiasing)
            config = gl.Config(sample_buffers=1, samples=4, \
                               depth_size=16, double_buffer=True)
            window = pyglet.window.Window(resizable=True, config=config)
        except pyglet.window.NoSuchConfigException:
            # Fall back to no multisampling for old hardware
            window = pyglet.window.Window(resizable=True)
        return window
    
        
    def glSetup(self):
        gl.glClearColor(0.2, 0.2, 0.3, 1)
        gl.glColor3f(1, 1, 1)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_CULL_FACE)
        gl.glEnable(gl.GL_NORMALIZE)
        
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, vec(50, 50, 10, 0))
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_SPECULAR, vec(5, 5, 10, 1))
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_DIFFUSE, vec(1, 1, 1, 1))
        gl.glLightfv(gl.GL_LIGHT1, gl.GL_AMBIENT, vec(0, 0, 0, 1.0))
        gl.glLightfv(gl.GL_LIGHT1, gl.GL_SPECULAR, vec(0.6, 0.6, 0.6, 1.0))
        gl.glLightfv(gl.GL_LIGHT1, gl.GL_DIFFUSE, vec(0.8, 0.8, 0.8, 1))
        gl.glLightfv(gl.GL_LIGHT1, gl.GL_POSITION, vec(-10.0, -20.0, 20.0, 0))

        # Create a Material and Group for the Model
        diffuse = [0.5, 0.5, 0.3, 1.0]
        ambient = [0.5, 0.5, 0.3, 1.0]
        specular = [1.0, 1.0, 1.0, 1.0]
        emission = [0.0, 0.0, 0.0, 1.0]
        shininess = 50
        material = pyglet.model.Material("", diffuse, ambient, specular, emission, shininess)
        self.group = pyglet.model.MaterialGroup(material=material)

        # Create a Material and Group for the ground plane
        diffuse = [0.02, 0.02, 0.023, 1.0]
        ambient = [0.01, 0.01, 0.01, 1.0]
        specular = [0.05, 0.05, 0.07, 1.0]
        emission = [0.0, 0.0, 0.0, 1.0]
        shininess = 10
        material2 = pyglet.model.Material("ground", diffuse, ambient, specular, emission, shininess)
        self.group2 = pyglet.model.MaterialGroup(material=material2)


    def makeBatches(self):
        # convert geometry lists to overall bounding box and
        # OpenGL-specific representations
        
        bbx = False
        # values for adding small offset to a point to make single-point bounding box
        epsP = point(epsilon,epsilon,epsilon)
        epsM = point(-epsilon,-epsilon,-epsilon)

        ## utility funtion to update bounding box based on a list of points
        def upbb(pp,box):
            for p in pp:
                if box:
                    box = [ [min(box[0][0],p[0]),min(box[0][1],p[1]),min(box[0][2],p[2])],
                            [max(box[1][0],p[0]),max(box[1][1],p[1]),max(box[1][2],p[2])] ]
                else:
                    box = [ add(p,epsM),add(p,epsP) ]
            return box
            
        if self.__points == [] and \
           self.__lines == [] and \
           self.__linestrips == [] and \
           self.__surfaces == []:
            raise ValueError('nothing to render')

        for p in self.__points:
            pp = p[0]
            pc = p[1]
            bbx = upbb(pp[1],bbx)
            self.__batch1.add(int(len(pp[1])/3),gl.GL_POINTS,self.group,pp,pc)

        for l in self.__lines:
            ll = l[0]
            lc = l[1]
            pp = [list(ll[1][0:3])]
            pp.append(list(ll[1][3:6]))
            bbx = upbb(pp,bbx)
            self.__batch1.add(int(len(ll[1])/3),gl.GL_LINES,self.group,ll,lc)
            
        for l in self.__linestrips:
            ll = l[0]
            lc = l[1]
            li = l[2]
            pp = []
            ln = int(len(ll[1])/3)
            for i in range(ln):
                pp.append(ll[1][i*3:(i+1)*3])
            bbx =upbb(pp,bbx)
            self.__batch1.add_indexed(int(len(ll[1])/3),
                                    gl.GL_LINES,
                                    self.group,
                                    li,ll,lc)

        for s in self.__surfaces:
            vert = s[0]
            norm = s[1]
            ind = s[2]
            pp = []
            ln = int(len(vert)/3)
            for i in range(ln):
                pp.append(vert[i*3:(i+1)*3])
            bbx = upbb(pp,bbx)
            self.__batch2.add_indexed(int(len(vert)/3),
                                    gl.GL_TRIANGLES,
                                    self.group,
                                    ind,
                                    ('v3f/static', vert),
                                    ('n3f/static', norm))

        ## Create a ground plane
        self.__batch3.add_indexed(4,
                                gl.GL_TRIANGLES,
                                self.group2,
                                [0,1,2,2,3,0],
                                ('v3f/static',[bbx[0][0]*1.2,bbx[0][1]*1.2,bbx[0][2]-0.1,
                                               bbx[1][0]*1.2,bbx[0][1]*1.2,bbx[0][2]-0.1,
                                               bbx[1][0]*1.2,bbx[1][1]*1.2,bbx[0][2]-0.1,
                                               bbx[0][0]*1.2,bbx[1][1]*1.2,bbx[0][2]-0.1]),
                                ('n3f/static',[0,0,1]*4))
                              
                                  
                                
        self.__bbox = bbx
        rnge = sub(self.__bbox[1],self.__bbox[0])
        mdim = max(rnge)
        mdim = max(mdim,10)

        # print("scene bounding box: ",vstr(bbx))

        
    def __init__(self):

        super().__init__()
        # these layers added for compatibility with ezdxf_drawable
        # default layer list.  For now, selecting a layer has
        # no effect on drawing
        self.layerlist = [False, '0','PATHS','DRILLS','DOCUMENTATION']

        def on_key_press(symbol,modifiers):
            if symbol == pyglet.window.key.UP:
                self.__cameradist *= 1.1
                if self.__cameradist > self.__maxcameradist:
                    self.__cameradist = self.__maxcameradist
            elif symbol == pyglet.window.key.DOWN:
                self.__cameradist *= .90909090
                if self.__cameradist < self.__mincameradist:
                    self.__cameradist = self.__mincameradist
            elif symbol == pyglet.window.key.RETURN:
                self.__cameradist = self.camerastartdist
                self.__rx = 0
                self.__ry = 0
            elif symbol == pyglet.window.key.P:
                self.__drawground = not self.__drawground
            elif symbol == pyglet.window.key.M:
                self.__legend = not self.__legend
            elif symbol == pyglet.window.key.L:
                if not (self.__light0 or self.__light1):
                    self.__light0 = True
                    gl.glEnable(gl.GL_LIGHT0)
                    gl.glDisable(gl.GL_LIGHT1)
                elif self.__light0 and not self.__light1:
                    self.__light1 = True
                    gl.glEnable(gl.GL_LIGHT0)
                    gl.glEnable(gl.GL_LIGHT1)
                elif self.__light0 and self.__light1:
                    self.__light0 = self.__light1 = False
            
                    
        self.__window = self.window()
        self.__window.push_handlers(on_key_press)
        self.__window.projection = pyglet.window.Projection3D(zfar=1000.0)
        self.glSetup()
        self.__center= point(0,0)
        self.__magnify = 0.08
        self.__arcres = 5
        self.__cameradist = self.camerastartdist = 100.0
        self.__maxcameradist = 900.0
        self.__mincameradist = 10.0
        self.__light0 = False
        self.__light1 = False
        self.__legend = True
        self.__drawground = True
        self.__rx = 0
        self.__ry = 0
        self.__points= []
        self.__lines= []
        self.__linestrips=[]
        self.__surfaces=[]
        self.__labels = []
        self.__batch1 = graphics.Batch() # use for lines, points, etc.
        self.__batch2 = graphics.Batch() # use for surfaces
        self.__batch3 = graphics.Batch() # use for environmental features, e.g. ground plane
        self.__bbox = False

        @self.__window.event
        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            self.__rx += dx
            self.__ry -= dy
            return pyglet.event.EVENT_HANDLED



        @self.__window.event
        def on_draw():
            self.__window.clear()
            gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
            
            self.__window.projection = pyglet.window.Projection3D(zfar=1000.0)
            # gl.glEnable(gl.GL_BLEND)
            gl.glColor3f(1., 1., 1.)

            #gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glLoadIdentity()
            cent = scale3(add(self.__bbox[0],self.__bbox[1]),0.5)
            rnge = sub(self.__bbox[1],self.__bbox[0])
            mdim = max(rnge)
            mdim = max(mdim,10)
            
            gl.glTranslatef(-cent[0],-cent[1],-cent[2]-1*self.__cameradist)
            #gl.glTranslatef(0,0,-1*self.__cameradist)
            gl.glRotatef(self.__ry%360.0, 1, 0, 0)
            gl.glRotatef(self.__rx%360.0, 0, 1, 0)

            #draw "ground plane"
            if self.__drawground:
                gl.glEnable(gl.GL_LIGHTING)
                gl.glEnable(gl.GL_LIGHT0)
                self.__batch3.draw()  

            gl.glDisable(gl.GL_LIGHTING)
            self.__batch1.draw()
            
            if self.__light0 or self.__light1:
                gl.glEnable(gl.GL_LIGHTING)
                self.__batch2.draw()

            #gl.glDisable(gl.GL_BLEND)
            gl.glDisable(gl.GL_LIGHTING)
            gl.glColor3f(1., 1., 1.)
            gl.glScalef(0.05, 0.05, 0.05)
            for label in self.__labels:
                gl.glPushMatrix()
                gl.glTranslatef(0,0,label[1])
                label[0].draw()
                gl.glPopMatrix()
                
            if self.__legend:
                gl.glMatrixMode(gl.GL_MODELVIEW)
                gl.glPushMatrix()
                gl.glLoadIdentity()
                gl.glMatrixMode(gl.GL_PROJECTION)
                gl.glPushMatrix()
                gl.glLoadIdentity()

                self.__window.projection = pyglet.window.Projection2D()
                label = pyglet.text.HTMLLabel(yapCAD_legend,
                                              x=10, y=30,
                                              anchor_x='left',anchor_y='bottom',
                                              width=400,
                                              multiline=True)
                label.draw()
                gl.glMatrixMode(gl.GL_PROJECTION)
                gl.glPopMatrix()
                gl.glMatrixMode(gl.GL_MODELVIEW)
                gl.glPopMatrix()
         
            return pyglet.event.EVENT_HANDLED
            

    def __repr__(self):
        return 'an instance of pygletDraw'

    ## properties

    @property
    def magnify(self):
        return self.__magnify

    def _set_magnify(self,mag):
        self.__magnify=mag

    @magnify.setter
    def magnify(self,mag=False):
        if not isinstance(mag,(int,float)):
            raise ValueError('invalid magnification ' + str(mag))
        if isinstance(mag,bool) and mag ==False:
            mag = 1.0
        elif mag < epsilon:
            mag = epsilon
        self._set_magnify(mag)

    @property
    def cameradist(self):
        return self.__cameradist

    def _set_cameradist(self,dist):
        self.__cameradist=dist

    @cameradist.setter
    def cameradist(self,dist=False):
        if not isinstance(dist,(int,float)):
            raise ValueError('invalid camera distance ' + str(dist))
        if isinstance(dist,bool) and dist ==False:
            dist = 100.0
        elif dist < self.__mincameradist:
            dist = self.__mincameradist
        self._set_cameradist(dist)

    
    ## Overload virtual yapcad.drawable base class dawing methods
    
    def draw_point(self,p):
        ar = self.__arcres
        self.__arcres=30
        super().draw_point(p)
        self.__arcres=ar

    def draw_line(self,p1,p2):
        color = self.thing2color(self.linecolor,'f')
        self.__lines.append([ ('v3f',(p1[0],p1[1],p1[2],
                                      p2[0],p2[1],p2[2])),
                              ('c3f',tuple(color + color)) ])
            
    def draw_arc(self,p,r,start,end):
        res = self.__arcres
        # resolution of arc sampling in degrees
        points = []
        if not (start==0 and end==360):
            start = start%360.0
            end = end % 360.0
            if end < start:
                end = end + 360
        theta = start*pi2/360.0
        points.append(add(p,[math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]))
        for a in range(round(start),round(end),res):
            theta = a*pi2/360.0
            pp = [math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]
            pp = add(pp,p)
            points.append(pp)
        theta = end*pi2/360.0
        points.append(add(p,[math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]))

        self.draw_linestrip(points)

    def draw_text(self,text,location,
                  align='left',
                  attr={'name': 'Times New Roman',
                        'size': 12}):
        name = 'Times New Roman'
        if 'font_name' in attr:
            name = attr['font_name']
        bold = False
        if 'bold' in attr:
            bold = attr['bold']
        italic = False
        if 'italic' in attr:
            italic = attr['italic']
        underline = False
        if 'underline' in attr:
            underline = attr['underline']
        size = 12
        if 'size' in attr:
            size = attr['size']
        elif 'height' in attr:
            size *= attr['height']
        anchor_x='left'
        if 'anchor_x' in attr:
            anchor_x = attr['anchor_x']
        anchor_y='center'
        if 'anchor_y' in attr:
            anchor_y = attr['anchor_y']
        col = []
        if 'color' in attr:
            col = attr['color']
        elif self.linecolor:
            col = self.thing2color(self.linecolor,'b')
            # print ("color: ",col)
        else:
            col = [255,255,255]
            
        color = tuple(col + [255])
            
        x = location[0]
        y = location[1]
        
        self.__labels.append([pyglet.text.Label(text,
                                                font_name=name,
                                                font_size=size*self.__magnify,
                                                x=x*20.0, #fudge factor to make font rendering look OK
                                                y=y*20.0,
                                                align=align.lower(),
                                                color=color,
                                                bold=bold,
                                                italic=italic,
                                                #underline=underline,
                                                anchor_x=anchor_x,
                                                anchor_y=anchor_y),
                              location[2]*20.0])
        
    ## OpenGL-specific drawing methods
    
    def draw_linestrip(self,points):
        # we simulate a linestrip useing GL_LINES and indexed drawing.
        # This prevents extra lines
        p2 = []
        i2 = []
        color = self.thing2color(self.linecolor,'f')
        for p in points:
            p2 = p2 + p[0:3]
        for i in range(1,len(points)):
            i2 = i2 + [i-1, i]
        c2 = color * int(len(p2)/3)
        self.__linestrips.append([ ('v3f', tuple(p2)),\
                                   ('c3f',tuple(c2)),
                                   tuple(i2) ])
            

    def draw_surface(self,points,normals,faces):
        vrts = []
        nrms= []
        inds = []
        for p in points:
            vrts+=p[0:3]
        for n in normals:
            nrms+=n[0:3]
        for f in faces:
            inds+=f[0:3]
        self.__surfaces.append([vrts,nrms,inds])    

    ## overload base-class virtual display method
    def display(self):
        self.makeBatches()
        pyglet.app.run()

