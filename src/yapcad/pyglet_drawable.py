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
from yapcad.geom3d import *
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

class Material():
    """
    Representation of an OpenGL-style material
    """

    def __repr__(self):
        return f"Material(ambient={self.ambient}, diffuse={self.diffuse}, specular={self.specular}, emission={self.emission}, shininess={self.shininess}, desc=\"{self.desc}\")"

    def __checkargs(self,v):
        if ( not (len(v) == 3 or len(v) == 4) or
             len(list(filter(lambda x: not isinstance(x,(int,float)),v))) > 0 or
             len(list(filter(lambda x: x < 0.0 or x > 1.0,v))) > 0):
            raise ValueError('bad arguments to property setter')
        if len(v) == 3:
            return [v[0],v[1],v[2],1.0]
        else:
            return v

    
    def __init__(self,**kwargs):
        self.__ambient = vec(0.2, 0.2, 0.2, 1.0)
        self.__diffuse = vec(0.5, 0.5, 0.5, 1.0)
        self.__specular = vec(0.5, 0.5, 0.5, 1.0)
        self.__emission = vec(0.0, 0.0, 0.0, 1.0)
        self.__shininess = 50
        self.__desc = ""
        if 'ambient' in kwargs:
            self.ambient = kwargs['ambient']
        if 'diffuse' in kwargs:
            self.diffuse = kwargs['diffuse']
        if 'specular' in kwargs:
            self.specular = kwargs['specular']
        if 'emission' in kwargs:
            self.emission = kwargs['emission']
        if 'shininess' in kwargs:
            self.shininess = kwargs['shininess']
        if 'desc' in kwargs:
            self.desc = kwargs['desc']

    @property
    def ambient(self):
        return list(self.__ambient)

    @ambient.setter
    def ambient(self,*v):
        if len(v) == 1:
            v = v[0]
        v = self.__checkargs(v)
        self.__ambient = vec(*v)
        
    @property
    def diffuse(self):
        return list(self.__diffuse)

    @diffuse.setter
    def diffuse(self,*v):
        if len(v) == 1:
            v = v[0]
        v = self.__checkargs(v)
        self.__diffuse = vec(*v)
        
    @property
    def specular(self):
        return list(self.__specular)

    @specular.setter
    def specular(self,*v):
        if len(v) == 1:
            v = v[0]
        v = self.__checkargs(v)
        self.__specular = vec(*v)

    @property
    def emission(self):
        return list(self.__emission)

    @emission.setter
    def emission(self,*v):
        if len(v) == 1:
            v = v[0]
        v = self.__checkargs(v)
        self.__emission = vec(*v)

    @property
    def shininess(self):
        return self.__shininess

    @shininess.setter
    def shininess(self,v):
        if not isinstance(v,(int,float)) or v < 0 or v > 128:
            raise ValueError('bad shininess')
        self.__shininess=v

    @property
    def desc(self):
        return self.__desc

    @desc.setter
    def desc(self,v):
        if not isinstance(v,str):
            raise ValueError('bad description')
        self.__desc = v

    @property
    def material(self):
        return pyglet.model.Material(self.desc,
                                     self.diffuse,
                                     self.ambient,
                                     self.specular,
                                     self.emission,
                                     self.shininess)

        

## global materials dictionary
materials = {}

materials['default'] = Material(ambient=[0.5, 0.3, 0.2, 1.0],
                                diffuse = [0.6, 0.5, 0.3, 1.0],
                                specular = [1.0, 0.8, 0.5, 1.0],
                                emission = [0.0, 0.0, 0.0, 1.0],
                                shininess = 20,
                                desc = "default material, dull gold")


materials['groundplane'] = Material(diffuse = [0.02, 0.02, 0.023, 1.0],
                                    ambient = [0.01, 0.01, 0.02, 1.0],
                                    specular = [0.05, 0.06, 0.08, 1.0],
                                    emission = [0.0, 0.0, 0.0, 1.0],
                                    shininess = 20,
                                    desc = "ground plane material")

## The following materials a from tables published at
## http://devernay.free.fr/cours/opengl/materials.html and
## http://www.it.hiof.no/~borres/j3d/explain/light/p-materials.html

materials['emerald'] = Material(ambient=[0.0215, 0.1745, 0.0215, 1.0],
                                diffuse=[0.07568, 0.61424, 0.07568, 1.0],
                                specular=[0.633, 0.727811, 0.633, 1.0],
                                emission = [0.0, 0.0, 0.0, 1.0],
                                shininess = 0.6*128,
                                desc= "emerald")

materials['jade'] = Material(ambient=[0.135, 0.2225, 0.1575, 1.0],
                             diffuse=[0.54, 0.89, 0.63, 1.0],
                             specular=[0.316228, 0.316228, 0.316228, 1.0],
                             emission = [0.0, 0.0, 0.0, 1.0],
                             shininess = 0.1*128)

materials['obsidian'] = Material(ambient=[0.05375, 0.05, 0.06625, 1.0],
                                 diffuse=[0.18275, 0.17, 0.22525, 1.0],
                                 specular=[0.332741, 0.328634, 0.346435, 1.0],
                                 emission = [0.0, 0.0, 0.0, 1.0],
                                 shininess = 0.3*128)

materials['pearl'] = Material(ambient=[0.25, 0.20725, 0.20725, 1.0],
                              diffuse=[0.9999, 0.829, 0.829, 1.0],
                              specular=[0.296648, 0.296648, 0.296648, 1.0],
                              emission = [0.0, 0.0, 0.0, 1.0],
                              shininess = 0.088*128)

materials['ruby'] = Material(ambient=[0.1745, 0.01175, 0.01175, 1.0],
                             diffuse=[0.61424, 0.04136, 0.04136, 1.0],
                             specular=[0.727811, 0.626959, 0.626959, 1.0],
                             emission = [0.0, 0.0, 0.0, 1.0],
                             shininess = 0.6*128)

materials['turquoise'] = Material(ambient=[0.1, 0.18725, 0.1745, 1.0],
                                  diffuse=[0.396, 0.74151, 0.69102, 1.0],
                                  specular=[0.297254, 0.30829, 0.306678, 1.0],
                                  emission = [0.0, 0.0, 0.0, 1.0],
                                  shininess = 0.1*128)

materials['brass'] = Material(ambient=[0.329412, 0.223529, 0.027451, 1.0],
                              diffuse=[0.780392, 0.568627, 0.113725, 1.0],
                              specular=[0.992157, 0.941176, 0.807843, 1.0],
                              emission = [0.0, 0.0, 0.0, 1.0],
                              shininess = 0.21794872*128)

materials['bronze'] = Material(ambient=[0.2125, 0.1275, 0.054, 1.0],
                               diffuse=[0.714, 0.4284, 0.18144, 1.0],
                               specular=[0.393548, 0.271906, 0.166721, 1.0],
                               emission = [0.0, 0.0, 0.0, 1.0],
                               shininess = 0.2*128)

materials['polished bronze'] = Material(
    ambient =[0.25, 0.148, 0.06475, 1.0],
    diffuse =[0.4, 0.2368, 0.1036, 1.0],
    specular =[0.774597, 0.458561, 0.200621, 1.0],
    shininess =76.8)

materials['chrome'] = Material(ambient=[0.25, 0.25, 0.25, 1.0],
                               diffuse=[0.4, 0.4, 0.4, 1.0],
                               specular=[0.774597, 0.774597, 0.774597, 1.0],
                               emission = [0.0, 0.0, 0.0, 1.0],
                               shininess = 0.6*128)

materials['copper'] = Material(ambient=[0.19125, 0.0735, 0.0225, 1.0],
                               diffuse=[0.7038, 0.27048, 0.0828, 1.0],
                               specular=[0.256777, 0.137622, 0.086014, 1.0],
                               emission = [0.0, 0.0, 0.0, 1.0],
                               shininess = 0.1*128)

materials['polished copper'] = Material(
    ambient =[ 0.2295, 0.08825, 0.0275, 1.0 ],
    diffuse =[0.5508, 0.2118, 0.066, 1.0 ],
    specular =[0.580594, 0.223257, 0.0695701, 1.0 ],
    shininess =51.2)
          
materials['gold'] =  Material( ambient=[0.24725, 0.1995, 0.0745, 1.0],
                               diffuse=[0.75164, 0.60648, 0.22648, 1.0],
                               specular=[0.628281, 0.555802, 0.366065, 1.0],
                               emission = [0.0, 0.0, 0.0, 1.0],
                               shininess = 0.4*128)

materials['polished gold'] = Material(
    ambient =[ 0.24725, 0.2245, 0.0645, 1.0 ],
    diffuse =[0.34615, 0.3143, 0.0903, 1.0 ],
    specular =[ 0.797357, 0.723991, 0.208006, 1.0],
    shininess =83.2)
          
materials['tin'] = Material(
    ambient =[ 0.105882, 0.058824, 0.113725, 1.0 ],
    diffuse =[0.427451, 0.470588, 0.541176, 1.0 ],
    specular =[0.333333, 0.333333, 0.521569, 1.0 ],
    shininess = 9.84615)

materials['silver'] = Material(
    ambient=[0.19225, 0.19225, 0.19225, 1.0],
    diffuse=[0.50754, 0.50754, 0.50754, 1.0],
    specular=[0.508273, 0.508273, 0.508273, 1.0],
    emission = [0.0, 0.0, 0.0, 1.0],
    shininess = 0.4*128)

materials['polished silver'] = Material(
    ambient =[ 0.23125, 0.23125, 0.23125, 1.0 ],
    diffuse =[0.2775, 0.2775, 0.2775, 1.0 ],
    specular =[0.773911, 0.773911, 0.773911, 1.0 ],
    shininess =89.6)

materials['black plastic'] = Material(ambient=[0.0, 0.0, 0.0, 1.0],
                                      diffuse=[0.01, 0.01, 0.01, 1.0],
                                      specular=[0.50, 0.50, 0.50, 1.0],
                                      shininess=.25*128)

materials['cyan plastic'] = Material(ambient=[0.0, 0.1, 0.06, 1.0],
                                     diffuse=[0.0, 0.50980392, 0.50980392, 1.0],
                                     specular=[0.50196078, 0.50196078,
                                               0.50196078, 1.0],
                                     shininess=.25*128)

materials['green plastic'] = Material(ambient=[0.0, 0.0, 0.0, 1.0],
                                      diffuse=[0.1, 0.35, 0.1, 1.0],
                                      specular=[0.45, 0.55, 0.45, 1.0],
                                      shininess=.25*128)

materials['red plastic'] = Material(ambient=[0.0, 0.0, 0.0, 1.0],
                                    diffuse=[0.5, 0.0, 0.0, 1.0],
                                    specular=[0.7, 0.6, 0.6, 1.0],
                                    emission=[0.0, 0.0, 0.0, 1.0],
                                    shininess=.25*128)

materials['white plastic'] = Material(ambient=[0.0, 0.0, 0.0, 1.0],
                                      diffuse=[0.55, 0.55, 0.55, 1.0],
                                      specular=[0.70, 0.70, 0.70, 1.0],
                                      shininess=.25*128)

materials['yellow plastic'] = Material(ambient=[0.0, 0.0, 0.0, 1.0],
                                       diffuse=[0.5, 0.5, 0.0, 1.0],
                                       specular=[0.60, 0.60, 0.50, 1.0],
                                       shininess=.25*128)

materials['black rubber'] = Material(ambient=[0.02, 0.02, 0.02, 1.0],
                                     diffuse=[0.01, 0.01, 0.01, 1.0],
                                     specular=[0.4, 0.4, 0.4, 1.0],
                                     emission=[0.0, 0.0, 0.0, 1.0],
                                     shininess=.078125*128)

materials['cyan rubber'] = Material(ambient=[0.0, 0.05, 0.05, 1.0],
                                    diffuse=[0.4, 0.5, 0.5, 1.0],
                                    specular=[0.04, 0.7, 0.7, 1.0],
                                    emission=[0.0, 0.0, 0.0, 1.0],
                                    shininess=.078125*128)

materials['green rubber'] = Material(ambient=[0.0, 0.05, 0.0, 1.0],
                                     diffuse=[0.4, 0.5, 0.4, 1.0],
                                     specular=[0.04, 0.7, 0.04, 1.0],
                                     shininess=.078125*128)

materials['red rubber'] = Material(ambient=[0.05, 0.0, 0.0, 1.0],
                                   diffuse=[0.5, 0.4, 0.4, 1.0],
                                   specular=[0.7, 0.04, 0.04, 1.0],
                                   shininess=.078125*128)

materials['white rubber'] = Material(ambient=[0.05, 0.05, 0.05, 1.0],
                                     diffuse=[0.5, 0.5, 0.5, 1.0],
                                     specular=[0.7, 0.7, 0.7, 1.0],
                                     emission=[0.0, 0.0, 0.0, 1.0],
                                     shininess=.078125*128)

materials['yellow rubber'] = Material(ambient=[0.05, 0.05, 0.0, 1.0],
                                      diffuse=[0.5, 0.5, 0.4, 1.0],
                                      specular=[0.7, 0.7, 0.04, 1.0],
                                      shininess=.078125*128)


## "throwaway" class for keeping track of object stuff
class GeomObject:
    pass

        
## class to provide openGL drawing functionality
class pygletDraw(drawable.Drawable):
    """
    yapCAD ``drawable`` subclass for OpenGL rendering with pyglet
    """

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
        self.__fps_display = pyglet.window.FPSDisplay(window)
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

        # Create a default material and Group for the Model
        material = materials['default'].material
        self.group = pyglet.model.MaterialGroup(material=material)

        # Create a Material and Group for the ground plane
        material2 = materials['groundplane'].material
        self.group2 = pyglet.model.MaterialGroup(material=material2)

    ## utility funtion to update bounding box based on a list of points
    def upbb(self,pp,box,offset=[0,0,0,1]):
        epsP = point(epsilon,epsilon,epsilon)
        epsM = point(-epsilon,-epsilon,-epsilon)
        for p in pp:
            p = add(p,offset)
            if box:
                box = [ [min(box[0][0],p[0]),min(box[0][1],p[1]),min(box[0][2],p[2])],
                        [max(box[1][0],p[0]),max(box[1][1],p[1]),max(box[1][2],p[2])] ]
            else:
                box = [ add(p,epsM),add(p,epsP) ]
        return box
        

    def addSurface(self,s,batch,group,offset,bbx):
        assert s[0] == 'surface'
        vert = s[1]
        norm = s[2]
        ind = s[3]
        pp = []
        ln = int(len(vert)/3)
        for i in range(ln):
            pp.append(vert[i*3:(i+1)*3])
        bbx = self.upbb(pp,bbx,offset)
        batch.add_indexed(int(len(vert)/3),
                          gl.GL_TRIANGLES,
                          group,
                          ind,
                          ('v3f/static', vert),
                          ('n3f/static', norm))
        return bbx

    __immediate = {}
        
    def register_immediate(self,func):
        """ register a function for immediate-mode drawing"""
        self.__immediate[func.__name__] = func

    def makeBatches(self):
        # convert geometry lists to overall bounding box and
        # OpenGL-specific representations
        
        bbx = False
        # values for adding small offset to a point to make single-point bounding box
        if ( self.__points == [] and 
             self.__lines == [] and 
             self.__linestrips == [] and 
             len(self.__objectdict) == 0 and 
             self.__surfaces == [] ):
            raise ValueError('nothing to render')

        for p in self.__points:
            pp = p[0]
            pc = p[1]
            bbx = self.upbb(pp[1],bbx)
            self.__batch1.add(int(len(pp[1])/3),gl.GL_POINTS,self.group,pp,pc)

        for l in self.__lines:
            ll = l[0]
            lc = l[1]
            pp = [list(ll[1][0:3])]
            pp.append(list(ll[1][3:6]))
            bbx = self.upbb(pp,bbx)
            self.__batch1.add(int(len(ll[1])/3),gl.GL_LINES,self.group,ll,lc)
            
        for l in self.__linestrips:
            ll = l[0]
            lc = l[1]
            li = l[2]
            pp = []
            ln = int(len(ll[1])/3)
            for i in range(ln):
                pp.append(ll[1][i*3:(i+1)*3])
            bbx =self.upbb(pp,bbx)
            self.__batch1.add_indexed(int(len(ll[1])/3),
                                    gl.GL_LINES,
                                    self.group,
                                    li,ll,lc)

            

        def lineSurface(s,obj):
            assert s[0] == 'surface'
            vert = s[1]
            ind = s[3]
            drawn = []
            for j in range(0,len(ind),3):
                i = ind[j:(j+3)]
                l1 = i[0],i[1]
                if l1[0] > l1[1]:
                    l1 = i[1],i[0]
                l2 = i[1],i[2]
                if l2[0] > l2[1]:
                    l2 = i[2],i[1]
                l3 = i[2],i[0]
                if l3[0] > l3[1]:
                    l3 = i[0],i[2]
                if not l1 in drawn:
                    self.draw_line(vert[l1[0]*3:l1[0]*3+3],
                                   vert[l1[1]*3:l1[1]*3+3],
                                   entity=obj)
                    drawn.append(l1)
                if not l2 in drawn:
                    self.draw_line(vert[l2[0]*3:l2[0]*3+3],
                                   vert[l2[1]*3:l2[1]*3+3],
                                   entity=obj)
                    drawn.append(l2)
                if not l3 in drawn:
                    self.draw_line(vert[l3[0]*3:l3[0]*3+3],
                                   vert[l3[1]*3:l3[1]*3+3],
                                   entity=obj)
                    drawn.append(l3)
            for l in obj.lines:
                ll = l[0]
                lc = l[1]
                obj.linebatch.add(2,gl.GL_LINES,self.group,ll,lc)

                    
        for o in self.__objectdict.values():
            for s in o.surfaces:
                bbx = self.addSurface(s,o.batch,o.group,
                                      point(o.x,o.y,o.z),bbx)
                lineSurface(s,o)
            for l in o.lines:
                ll = l[0]
                lc = l[1]
                o.linebatch.add(2,gl.GL_LINES,self.group,ll,lc)

                

        for s in self.__surfaces:
            bbx = self.addSurface(s,self.__batch2,self.group,
                                  point(0,0,0),bbx)

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
            handled = False
            if symbol == pyglet.window.key.UP:
                handled = True
                self.__cameradist *= 1.1
                if self.__cameradist > self.__maxcameradist:
                    self.__cameradist = self.__maxcameradist
            elif symbol == pyglet.window.key.DOWN:
                handled = True
                self.__cameradist *= .90909090
                if self.__cameradist < self.__mincameradist:
                    self.__cameradist = self.__mincameradist
            elif symbol == pyglet.window.key.RETURN:
                handled = True
                self.__cameradist = self.camerastartdist
                self.__rx = 0
                self.__ry = 0
            elif symbol == pyglet.window.key.P:
                handled = True
                self.__drawground = not self.__drawground
            elif symbol == pyglet.window.key.M:
                handled = True
                self.__legend = not self.__legend
            elif symbol == pyglet.window.key.L:
                handled = True
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
            if handled:
                return pyglet.event.EVENT_HANDLED
            else:
                return
            
                    
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
        self.__light0 = True
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
        self.__batch2 = graphics.Batch() # default, use for surfaces
        self.__batch3 = graphics.Batch() # use for environmental features, e.g. ground plane
        self.__objectdict = {} # use for animated objects
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
            if not self.__bbox:
                self.__bbox =[point(-1,-1,-1),point(1,1,1)]
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

            # execute any immediate-mode registered rendering functions
            for f in self.__immediate.values():
                f()
                
            for obj in self.__objectdict.values():
                gl.glMatrixMode(gl.GL_MODELVIEW)
                gl.glPushMatrix()
                gl.glTranslatef(obj.x,obj.y,obj.z)
                gl.glRotatef(obj.rx, 1, 0, 0)
                gl.glRotatef(obj.ry, 0, 1, 0)
                gl.glRotatef(obj.rz, 0, 0, 1)
                if obj.lighting and (self.__light0 or self.__light1):
                    gl.glEnable(gl.GL_LIGHTING)
                    obj.batch.draw()
                else:
                    gl.glDisable(gl.GL_LIGHTING)
                    obj.linebatch.draw()
                gl.glPopMatrix()

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
                self.__fps_display.draw()
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


    @property
    def rx(self):
        return self.__rx

    @rx.setter
    def rx(self,v):
        if not isinstance(v,(int,float)):
            raise ValueError('invalid camera rotation')
        self.__rx = v%360.0

    @property
    def ry(self):
        return self.__ry

    @ry.setter
    def ry(self,v):
        if not isinstance(v,(int,float)):
            raise ValueError('invalid camera rotation')
        self.__ry = v%360.0

    @property
    def objectdict(self):
        return self.__objectdict

    
    ## Overload virtual yapcad.drawable base class dawing methods
    
    def draw_point(self,p):
        ar = self.__arcres
        self.__arcres=30
        super().draw_point(p)
        self.__arcres=ar

    def draw_line(self,p1,p2,entity=None,c1=None,c2=None):

        if not entity:
            entity=self
            elist = entity.__lines
        else:
            elist = entity.lines
        color1 = []
        color2 = []
        if not c1:
            color1 = self.thing2color(entity.linecolor,'f')
        else:
            color1 = self.thing2color(c1,'f')
        if not c2:
            color2 = self.thing2color(entity.linecolor,'f')
        else:
            color2 = self.thing2color(c2,'f')
        elist.append([ ('v3f',(p1[0],p1[1],p1[2],
                               p2[0],p2[1],p2[2])),
                       ('c3f',tuple(color1 + color2)) ])
            
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

    def make_object(self,name,**kwargs):
        """
        Create a new three-dmensional object with specified
        material properties and 6DOF in world-space that can be
        animated
        """
        if name in self.__objectdict.keys():
            raise ValueError('duplicate object key')
        obj = GeomObject()

        if 'lighting' in kwargs:
            if kwargs['lighting']:
                obj.lighting = True
            else:
                obj.lighting=False
        if 'position' in kwargs:
            x = kwargs['position']
            if not ispoint(x):
                raise ValueError('bad position')
            obj.x = x[0]
            obj.y = x[1]
            obj.z = x[2]
        else:
            obj.x = 0
            obj.y = 0
            obj.z = 0
        if 'rotation' in kwargs:
            r = kwargs['rotation']
            obj.rx = r[0]
            obj.ry = r[1]
            obj.rz = r[2]
        else:
            obj.rx = 0
            obj.ry = 0
            obj.rz = 0
        if 'material' in kwargs:
            mat = kwargs['material']
            if mat in materials:
                obj.material = materials[mat]
                m = obj.material.material
                obj.group = pyglet.model.MaterialGroup(material=m)
            else:
                raise ValueError('unknown material')
        else:
            obj.material=materials['default']
            obj.group = self.group

        if 'animate' in kwargs:
            pyglet.clock.schedule_interval(kwargs['animate'],1/60)

        if 'linecolor' in kwargs:
            obj.linecolor = kwargs['linecolor']
        else:
            obj.linecolor = 'gray'

        obj.batch = graphics.Batch()
        obj.linebatch = graphics.Batch()
        obj.surfaces = []
        obj.lines=[]
        

        self.__objectdict[name] = obj
        
    def draw_surface(self,points,normals=None,faces=None,name=None):
        """
        render a surface, either as part of default scene geometry
        (if ``name`` is not specified) or as a part of a named object.
        """
        if not (normals or faces):
            if not issurface(points):
                raise ValueError('bad surface passed to draw_surface')
            normals = points[2]
            faces = points[3]
            points = points[1]
        vrts = []
        nrms= []
        inds = []
        for p in points:
            vrts+=p[0:3]
        for n in normals:
            nrms+=n[0:3]
        for f in faces:
            inds+=f[0:3]
        if not name:
            self.__surfaces.append(['surface',vrts,nrms,inds])
        else:
            obj = self.__objectdict[name]
            obj.surfaces.append(['surface',vrts,nrms,inds])

    def draw_solid(self,solid,name=None):
        """
        Render a solid, either as part of the secene geometry
        """
        if not issolid(solid):
            raise ValueError('bad solid passed to draw_solid')
        for surface in solid[1]:
            self.draw_surface(surface,name=name)

    ## override base-class draw method
    def draw(self,x,name=None):
        if issolid(x):
            self.draw_solid(x,name=name)
        elif issurface(x):
            self.draw_surface(x,name=name)
        else:
            super().draw(x)

    ## override base-class virtual display method
    def display(self):
        self.makeBatches()
        pyglet.app.run()

