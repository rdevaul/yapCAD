## simple object-oriented framework for dxf-rendered drawable objects
## in yapCAD

import math
from geom import *
import pyglet
import pyglet.gl as gl
import pyglet.graphics as graphics
import drawable

## base class to provide dxf drawing functionality
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
    
    
    def camera(self):
        o = self.center
        gl.glOrtho(o[0] - 1. / self.magnify,
                   o[0] + 1. / self.magnify,
                   o[1] - 1. / self.magnify,
                   o[1] + 1. / self.magnify, -1., 1.)

    def glSetup(self):
        gl.glClearColor(0, 0, 0, 1)
        gl.glColor3f(1, 1, 1)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_CULL_FACE)
        #gl.glEnable(gl.GL_LIGHTING)
        #gl.glEnable(gl.GL_LIGHT0)
        #gl.glEnable(gl.GL_LIGHT1)

        # Create a Material and Group for the Model
        diffuse = [0.5, 0.5, 0.3, 1.0]
        ambient = [0.5, 0.5, 0.3, 1.0]
        specular = [1.0, 1.0, 1.0, 1.0]
        emission = [0.0, 0.0, 0.0, 1.0]
        shininess = 50
        material = pyglet.model.Material("", diffuse, ambient, specular, emission, shininess)
        self.group = pyglet.model.MaterialGroup(material=material)

    
    def __init__(self):
        self.window = self.window()
        self.window.projection = pyglet.window.Projection3D()
        self.glSetup()
        self.center= point(0,0)
        self.magnify = 0.08
        self.points= []
        self.rx = 0
        self.ry = 0
        self.lines= []
        self.linestrips=[]
        self.labels = []
        self.linecolor = 'white'

        @self.window.event
        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            self.rx += dx
            self.ry -= dy

        
        @self.window.event
        def on_draw():
            self.window.clear()
            gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
            #gl.glTranslatef(0, 0, -4)
            
            # gl.glEnable(gl.GL_BLEND)
            gl.glColor3f(1., 1., 1.)
            batch = graphics.Batch()
            for p in self.points:
                pp = p[0]
                pc = p[1]
                batch.add(int(len(pp[1])/3),gl.GL_POINTS,self.group,pp,pc)
            for l in self.lines:
                ll = l[0]
                lc = l[1]
                batch.add(int(len(ll[1])/3),gl.GL_LINES,self.group,ll,lc)
            for l in self.linestrips:
                ll = l[0]
                lc = l[1]
                li = l[2]
                batch.add_indexed(int(len(ll[1])/3),
                                  gl.GL_LINES,
                                  self.group,
                                  li,ll,lc)
            #gl.glPopMatrix()

            if False:
                gl.glMatrixMode(gl.GL_PROJECTION)
                gl.glLoadIdentity()
                self.camera()

            #gl.glMatrixMode(gl.GL_MODELVIEW)
            gl.glLoadIdentity()
            gl.glTranslatef(0,0,-100.)
            gl.glRotatef(self.ry%360.0, 1, 0, 0)
            gl.glRotatef(self.rx%360.0, 0, 1, 0)
            batch.draw()
            #gl.glDisable(gl.GL_BLEND)
            for label in self.labels:
                label.draw()
            
        super().__init__()

    def __repr__(self):
        return 'an instance of pygletDraw'

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
        res=5                   # resolution of arc sampling in degrees
        points = []
        if not (start==0 and end==360):
            start = start%360.0
            end = end % 360.0
            if end < start:
                end = end + 360
        
        for a in range(round(start),round(end),res):
            theta = a*pi2/360.0
            pp = [math.cos(theta)*r,math.sin(theta)*r,0.0,1.0]
            pp = add(pp,p)
            points.append(pp)

        self.draw_linestrip(points)
        #self.lines.append([('v3f',tuple(verts)),
        #                   ('c3f',tuple(colors))])


    def draw_text(self,text,location,
                  attr={'name': 'Times New Roman',
                        'size': 12}):
        name = 'Times New Roman'
        if 'name' in attr:
            name = attr['name']
        size = 12
        if 'size' in attr:
            size = attr['size']
        anchor_x='center'
        if 'anchor_x' in attr:
            anchor_x = attr['anchor_x']
        anchor_y='center'
        if 'anchor_y' in attr:
            anchor_y = attr['anchor_y']
        x = location[0]
        y = location[1]
        
        self.labels.append(pyglet.text.Label(text,
                                             font_name=name,
                                             font_size=size*self.magnify,
                                             x=x,
                                             y=y,
                                             anchor_x=anchor_x,
                                             anchor_y=anchor_y))
        
    def display(self):
        pyglet.app.run()

