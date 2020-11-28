## simple yapCAD framework for openGL drawing using pyglet
## package
## Original Author: Richard W. DeVaul

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
            
        if self.points == [] and \
           self.lines == [] and \
           self.linestrips == [] and \
           self.surfaces == []:
            raise ValueError('nothing to render')

        for p in self.points:
            pp = p[0]
            pc = p[1]
            bbx = upbb(pp[1],bbx)
            self.batch1.add(int(len(pp[1])/3),gl.GL_POINTS,self.group,pp,pc)

        for l in self.lines:
            ll = l[0]
            lc = l[1]
            pp = [list(ll[1][0:3])]
            pp.append(list(ll[1][3:6]))
            bbx = upbb(pp,bbx)
            self.batch1.add(int(len(ll[1])/3),gl.GL_LINES,self.group,ll,lc)
            
        for l in self.linestrips:
            ll = l[0]
            lc = l[1]
            li = l[2]
            pp = []
            ln = int(len(ll[1])/3)
            for i in range(ln):
                pp.append(ll[1][i*3:(i+1)*3])
            bbx =upbb(pp,bbx)
            self.batch1.add_indexed(int(len(ll[1])/3),
                                    gl.GL_LINES,
                                    self.group,
                                    li,ll,lc)

        for s in self.surfaces:
            vert = s[0]
            norm = s[1]
            ind = s[2]
            pp = []
            ln = int(len(vert)/3)
            for i in range(ln):
                pp.append(vert[i*3:(i+1)*3])
            bbx = upbb(pp,bbx)
            self.batch2.add_indexed(int(len(vert)/3),
                                    gl.GL_TRIANGLES,
                                    self.group,
                                    ind,
                                    ('v3f/static', vert),
                                    ('n3f/static', norm))

        ## Create a ground plane
        self.batch3.add_indexed(4,
                                gl.GL_TRIANGLES,
                                self.group2,
                                [0,1,2,2,3,0],
                                ('v3f/static',[bbx[0][0]*1.2,bbx[0][1]*1.2,bbx[0][2]-0.1,
                                               bbx[1][0]*1.2,bbx[0][1]*1.2,bbx[0][2]-0.1,
                                               bbx[1][0]*1.2,bbx[1][1]*1.2,bbx[0][2]-0.1,
                                               bbx[0][0]*1.2,bbx[1][1]*1.2,bbx[0][2]-0.1]),
                                ('n3f/static',[0,0,1]*4))
                              
                                  
                                
        self.bbox = bbx
        rnge = sub(self.bbox[1],self.bbox[0])
        mdim = max(rnge)
        mdim = max(mdim,10)

        # print("scene bounding box: ",vstr(bbx))

        
    def __init__(self):

        def on_key_press(symbol,modifiers):
            if symbol == pyglet.window.key.UP:
                self.cameradist *= 1.1
                if self.cameradist > self.maxcameradist:
                    self.cameradist = self.maxcameradist
            elif symbol == pyglet.window.key.DOWN:
                self.cameradist *= .90909090
                if self.cameradist < self.mincameradist:
                    self.cameradist = self.mincameradist
            elif symbol == pyglet.window.key.RETURN:
                self.cameradist = self.camerastartdist
                self.rx = 0
                self.ry = 0
            elif symbol == pyglet.window.key.P:
                self.drawground = not self.drawground
            elif symbol == pyglet.window.key.M:
                self.legend = not self.legend
            elif symbol == pyglet.window.key.L:
                if not (self.light0 or self.light1):
                    self.light0 = True
                    gl.glEnable(gl.GL_LIGHT0)
                    gl.glDisable(gl.GL_LIGHT1)
                elif self.light0 and not self.light1:
                    self.light1 = True
                    gl.glEnable(gl.GL_LIGHT0)
                    gl.glEnable(gl.GL_LIGHT1)
                elif self.light0 and self.light1:
                    self.light0 = self.light1 = False
            
                    
        self.window = self.window()
        self.window.push_handlers(on_key_press)
        self.window.projection = pyglet.window.Projection3D()
        self.glSetup()
        self.center= point(0,0)
        self.magnify = 0.08
        self.arcres = 5
        self.cameradist = self.camerastartdist = 100.0
        self.maxcameradist = 255.0
        self.mincameradist = 10.0
        self.light0 = False
        self.light1 = False
        self.legend = True
        self.drawground = True
        self.points= []
        self.rx = 0
        self.ry = 0
        self.lines= []
        self.linestrips=[]
        self.surfaces=[]
        self.labels = []
        self.batch1 = graphics.Batch() # use for lines, points, etc.
        self.batch2 = graphics.Batch() # use for surfaces
        self.batch3 = graphics.Batch() # use for environmental features, e.g. ground plane
        self.linecolor = 'white'
        self.bbox = False

        @self.window.event
        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            self.rx += dx
            self.ry -= dy
            return pyglet.event.EVENT_HANDLED



        @self.window.event
        def on_draw():
            self.window.clear()
            gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
            
            self.window.projection = pyglet.window.Projection3D()
            # gl.glEnable(gl.GL_BLEND)
            gl.glColor3f(1., 1., 1.)

            #gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glLoadIdentity()
            cent = scale3(add(self.bbox[0],self.bbox[1]),0.5)
            rnge = sub(self.bbox[1],self.bbox[0])
            mdim = max(rnge)
            mdim = max(mdim,10)
            
            gl.glTranslatef(-cent[0],-cent[1],-cent[2]-1*self.cameradist)
            #gl.glTranslatef(0,0,-1*self.cameradist)
            gl.glRotatef(self.ry%360.0, 1, 0, 0)
            gl.glRotatef(self.rx%360.0, 0, 1, 0)

            #draw "ground plane"
            if self.drawground:
                gl.glEnable(gl.GL_LIGHTING)
                gl.glEnable(gl.GL_LIGHT0)
                self.batch3.draw()  

            gl.glDisable(gl.GL_LIGHTING)
            self.batch1.draw()
            
            if self.light0 or self.light1:
                gl.glEnable(gl.GL_LIGHTING)
                self.batch2.draw()

            #gl.glDisable(gl.GL_BLEND)
            gl.glDisable(gl.GL_LIGHTING)
            gl.glColor3f(1., 1., 1.)
            gl.glScalef(0.05, 0.05, 0.05)
            for label in self.labels:
                label.draw()
                
            if self.legend:
                gl.glMatrixMode(gl.GL_MODELVIEW)
                gl.glPushMatrix()
                gl.glLoadIdentity()
                gl.glMatrixMode(gl.GL_PROJECTION)
                gl.glPushMatrix()
                gl.glLoadIdentity()

                self.window.projection = pyglet.window.Projection2D()
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
            
        super().__init__()

    def __repr__(self):
        return 'an instance of pygletDraw'

    def draw_point(self,p):
        ar = self.arcres
        self.arcres=30
        super().draw_point(p)
        self.arcres=ar

    def draw_line(self,p1,p2):
        color = self.thing2color(self.linecolor,'f')
        self.lines.append([ ('v3f',(p1[0],p1[1],p1[2],
                                    p2[0],p2[1],p2[2])),
                            ('c3f',tuple(color + color)) ])
            
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
        self.linestrips.append([ ('v3f', tuple(p2)),\
                                 ('c3f',tuple(c2)),
                                 tuple(i2) ])
            

    def draw_arc(self,p,r,start,end):
        res = self.arcres
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
        #self.lines.append([('v3f',tuple(verts)),
        #                   ('c3f',tuple(colors))])

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
        #print ("len(points): ",len(points),
        #       "len(normals): ",len(normals),
        #       "len(faces): ",len(faces))
        self.surfaces.append([vrts,nrms,inds])


    
    def draw_text(self,text,location,
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
            print ("color: ",col)
        else:
            col = [255,255,255]
            
        color = tuple(col + [255])
            
        x = location[0]
        y = location[1]
        
        self.labels.append(pyglet.text.Label(text,
                                             font_name=name,
                                             font_size=size*self.magnify,
                                             x=x*20.0, #fudge factor to make font rendering look OK
                                             y=y*20.0,
                                             align='left',
                                             color=color,
                                             bold=bold,
                                             italic=italic,
                                             #underline=underline,
                                             anchor_x=anchor_x,
                                             anchor_y=anchor_y))
        
    def display(self):
        self.makeBatches()
        pyglet.app.run()

