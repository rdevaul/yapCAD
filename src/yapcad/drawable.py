## base class of drawable for yapCAD
## Copyright (c) 2020 Richard W. DeVaul
## Copyright (c) 2020 yapCAD contributors
## All rights reserved
## See licensing terms here: https://github.com/rdevaul/yapCAD/blob/master/LICENSE

from yapcad.geom import *

## Generic drawing functions -- assumed to use current coordinate
## transform and drawing pen (color, line weight, etc.)
 
class Drawable:
    """Base class for yapCAD drawables"""

    ## utility functions for drawing
    ## -----------------------------
    
    ## pure virtual functions -- override for specific rendering
    ## system
    def draw_arc(self,p,r,start,end):
        print("pure virtual draw_arc called: {}, {}, {}, {}".format(p,r,start,end)) 
        return

    def draw_line(self,p1,p2):
        print("pure virtual draw_line called: {}, {}".format(p1,p2))
        return

    def draw_text(self,text,location,
                  align='left',
                  attr={}):
        print("pure virtual draw_text called: {}, {}, {}, {}".format(text,location,align,attr))
        return

    ## non-virtual utility drawing functions 
    def draw_circle(self,p,r):
        self.draw_arc(p,r,0.0,360.0)
        return

    ## utility function to draw an "x", centered on point p inside square of
    ## dimension d
    def draw_x(self,p,d):
        hd=d/2
        p1=add(p,vect([-hd,-hd,0]))
        p2=add(p,vect([hd,hd,0]))
        p3=add(p,vect([-hd,hd,0]))
        p4=add(p,vect([hd,-hd,0]))
        self.draw_line(p1,p2)
        self.draw_line(p3,p4)
        return;

    ## utiitly function to draw a point based on the current point style
    def draw_point(self,p):
        if 'x' in self.__pointstyle:
            self.draw_x(p,self.__pointsize*2)

        if 'o' in self.__pointstyle:
            self.draw_circle(p,self.__pointsize)

    def __init__(self):
        self.__pointstyle = 'xo'
        self.__pointsize = 0.1
        self.__linetype = False
        self.__linewidth = 0.1
        self.__linecolor = False
        self.__fillcolor = False
        self.__polystyle = 'lines'
        self.__layer = False
        self.__layerlist = [ False, 'default' ]


    ## Various property functions

    @property
    def layerlist(self):
        return self.__layerlist

    def _set_layerlist(self,lst):
        self.__layerlist = lst

    @layerlist.setter
    def layerlist(self,lst):
        if isinstance(lst,list):
            self._set_layerlist(lst)
        else:
            raise ValueError('bad layer list ' + str(lst))
            
    
    @property
    def layer(self):
        return self.__layer
    
    def _set_layer(self,lyr):
        self.__layer = lyr
        
    @layer.setter
    def layer(self,lyr=False):
        if lyr in self.layerlist:
            self._set_layer(lyr)
        else:
            raise ValueError('bad layer: ' + str(lyr))
        
    @property
    def polystyle(self):
        return self.__polystyle

    def _set_polystyle(self,pst):
        self.__polystyle = pst
        
    @polystyle.setter
    def polystyle(self,pst=False):
        if pst in [False, 'points','lines','both']:
            self._set_polystyle(pst)
        else:
            raise ValueError('bad polystyle')
        
    @property
    def pointstyle(self):
        return self.__pointstyle

    def _set_pointstyle(self,pst):
        self.__pointstyle=pst
    
    @pointstyle.setter
    def pointstyle(self,pst=False):
        if pst in [False, 'x','o','xo']:
            self._set_pointstyle(pst)
        else:
            raise ValueError('bad pointstyle: ' + str(pst))

    @property
    def pointsize(self):
        return self.__pointsize

    def _set_pointsize(self,ps):
        self.__pointsize=ps
        
    @pointsize.setter
    def pointsize(self,ps=False):
        if not isinstance(ps,(int,float)):
            raise ValueError('invalid pointsize ' + str(ps))
        if isinstance(ps,bool) and ps == False:
            self._set_pointsize(0.1)
        elif ps < epsilon:
            ps = epsilon
        self._set_pointsize(ps)

    @property
    def linewidth(self):
        return self.__linewidth

    def _set_linewidth(self,lw):
        self.__linewidth=lw

    @linewidth.setter
    def linewidth(self,lw=False):
        if not isinstance(lw,(int,float)):
            raise ValueError('invalid pointsize ' + str(lw))
        if isinstance(lw,bool) and lw ==False:
            lw = 0.1
        elif lw < epsilon:
            lw = epsilon
        self._set_linewidth(lw)

    @property
    def linetype(self):
        return self.__linetype

    def _set_linetype(self,lt):
        self.__linetype = lt

    @linetype.setter
    def linetype(self,lt=False):
        if lt in [ False, 'Continuous' ]:       # only continuous linetype supported in base class
            self._set_linetype(False)
        else:
            raise ValueError('unsupported linetype ' + str(lt))

    ## color is a complex property -- it can be set as an AUTOCAD
    ## index color, a standard color name, or an RGB tripple

    def __checkcolor(self,c):
        def isbyte(x):
            return isinstance(x,int) and x >= 0 and x <= 255
        if isinstance(c,bool) and c == False:
            return True
        elif (isinstance(c,list) and len(c) == 3 \
              and isbyte(c[0]) and isbyte(c[1]) and isbyte(c[2])) or \
             (isinstance(c,int) and c >= 0 and c < len(self.colormapAUTOCAD)) or\
             isinstance(c,str) and c in self.colordict:
            return True
        return False
    
    ## line color
    @property
    def linecolor(self):
        return self.__linecolor

    def _set_linecolor(self,c):
        self.__linecolor=c

    @linecolor.setter
    def linecolor(self,c=False):
        if self.__checkcolor(c):
            self._set_linecolor(c)
        else:
            raise ValueError('bad linecolor ' + str(c))

    ## fill color
    @property
    def fillcolor(self):
        return self.__fillcolor

    def _set_fillcolor(self,c):
        self.__fillcolor = c

    @fillcolor.setter
    def fillcolor(self,c=False):
        if self.__checkcolor(c):
            self._set_fillcolor(c)
        else:
            raise ValueError('bad fillcolor ' + str(c))

    ## non-property methods
    
    def __repr__(self):
        return 'an abstract Drawable instance'

    def str(self):
        return self.__repr__()
    
    def print(self):
        print(self.str())

    def draw(self,x):
        if ispoint(x):
            self.draw_point(x)
        elif isline(x):
            self.draw_line(x[0],x[1])
        elif isarc(x):
            self.draw_arc(x[0],x[1][0],x[1][1],x[1][2])
        elif ispoly(x):
            if self.polystyle in ('points','lines','both'):
                if self.polystyle == 'points' or \
                   self.polystyle == 'both':
                    for e in x:
                        self.draw(e)
                if self.polystyle == 'lines' or \
                   self.polystyle == 'both':
                    for i in range(1,len(x)):
                        self.draw(line(x[i-1],x[i]))
            else:
                raise ValueError("bad value for polystyle: {}".format(self.polystyle))
            
        elif isgeomlist(x):
            for e in x:
                self.draw(e)
        else:
            raise ValueError('bad argument to Drawable.draw(): '.format(x))
        
    ## cause drawing page to be rendered -- pure virtual in base class
    def display(self):
        print('pure virtual display function called')
        return True


    ## Standard color names, ala HTML: https://www.rapidtables.com/web/color/RGB_Color.html
    colordict= {
        'black':   [ [0,0,0], 0],
        'white':   [ [255,255,255], 7],
        'red':     [ [255,0,0], 1],
        'lime':    [ [0,255,0],3],
        'blue':    [ [0,0,255], 5],
        'yellow':  [ [255,255,0], 2],
        'silver':  [ [192,192,192], 9],
        'gray':    [ [128,128,128], 8],
        'maroon':  [ [127,0,0], 14],
        'olive':   [ [127,127,0], 54],
        'green':   [ [0,127,0], 94],
        'aqua':    [ [0,255,255], 4],
        'teal':    [ [0,127,127], 134],
        'navy':    [ [0,0,128], 174],
        'fuchsia': [ [255,0,255], 6],
        'purple':  [ [128,0,128], 214] }
    
    ## AUTOCAD color index:
    ## https://forums.autodesk.com/t5/visual-lisp-autolisp-and-general/
    ## rgb-values-of-the-256-standard-indexed-colors/td-p/848912
    colormapAUTOCAD = [
        [0,0,0],
	[255,0,0],
	[255,255,0],
	[0,255,0],
	[0,255,255],
	[0,0,255],
	[255,0,255],
	[255,255,255],
	[128,128,128],
	[192,192,192],

	[255,0,0],
	[255,127,127],
	[165,0,0],
	[165,82,82],
	[127,0,0],
	[127,63,63],
	[76,0,0],
	[76,38,38],
	[38,0,0],
	[38,19,19],

	[255,63,0],
	[255,159,127],
	[165,41,0],
	[165,103,82],
	[127,31,0],
	[127,79,63],
	[76,19,0],
	[76,47,38],
	[38,9,0],
	[38,23,19],

	[255,127,0],
	[255,191,127],
	[168,82,0],
	[165,124,82],
	[127,63,0],
	[127,95,63],
	[76,38,0],
	[76,57,38],
	[38,19,0],
	[38,28,19],

	[255,191,0],
	[255,223,127],
	[165,124,0],
	[165,145,82],
	[127,95,0],
	[127,111,63],
	[76,57,0],
	[76,66,38],
	[38,28,0],
	[38,33,19],

	[255,255,0],
	[255,255,127],
	[165,165,0],
	[165,165,82],
	[127,127,0],
	[127,127,63],
	[76,76,0],
	[76,76,38],
	[38,38,0],
	[38,38,19],

	[191,255,0],
	[223,255,127],
	[124,165,0],
	[145,165,82],
	[95,127,0],
	[111,127,63],
	[57,76,0],
	[66,76,38],
	[28,38,0],
	[33,38,19],

	[127,255,0],
	[191,255,127],
	[82,165,0],
	[124,165,82],
	[63,127,0],
	[95,127,63],
	[38,76,0],
	[57,76,38],
	[19,38,0],
	[28,38,19],

	[63,255,0],
	[159,255,127],
	[41,165,0],
	[103,165,82],
	[31,127,0],
	[79,127,63],
	[19,76,0],
	[47,76,38],
	[9,38,0],
	[23,38,19],

	[0,255,0],
	[127,255,127],
	[0,165,0],
	[82,165,82],
	[0,127,0],
	[63,127,63],
	[0,76,0],
	[38,76,38],
	[0,38,0],
	[19,38,19],

	[0,255,63],
	[127,255,159],
	[0,165,41],
	[82,165,103],
	[0,127,31],
	[63,127,79],
	[0,76,19],
	[38,76,47],
	[0,38,9],
	[19,38,23],

	[0,255,127],
	[127,255,191],
	[0,165,82],
	[82,165,124],
	[0,127,63],
	[63,127,95],
	[0,76,38],
	[38,76,57],
	[0,38,19],
	[19,38,28],

	[0,255,191],
	[127,255,223],
	[0,165,124],
	[82,165,145],
	[0,127,95],
	[63,127,111],
	[0,76,57],
	[38,76,66],
	[0,38,28],
	[19,38,33],

	[0,255,255],
	[127,255,255],
	[0,165,165],
	[82,165,165],
	[0,127,127],
	[63,127,127],
	[0,76,76],
	[38,76,76],
	[0,38,38],
	[19,38,38],

	[0,191,255],
	[127,223,255],
	[0,124,165],
	[82,145,165],
	[0,95,127],
	[63,111,127],
	[0,57,76],
 	[38,66,76],
	[0,28,38],
	[19,33,38],

	[0,127,255],
	[127,191,255],
	[0,82,165],
	[82,124,165],
	[0,63,127],
	[63,95,127],
	[0,38,76],
	[38,57,76],
	[0,19,38],
	[19,28,38],

	[0,63,255],
	[127,159,255],
	[0,41,165],
	[82,103,165],
	[0,31,127],
	[63,79,127],
	[0,19,76],
	[38,47,76],
	[0,9,38],
	[19,23,38],

	[0,0,255],
	[127,127,255],
	[0,0,165],
	[82,82,165],
	[0,0,127],
	[63,63,127],
	[0,0,76],
	[38,38,76],
	[0,0,38],
	[19,19,38],

	[63,0,255],
	[159,127,255],
	[41,0,165],
	[103,82,165],
	[31,0,127],
	[79,63,127],
	[19,0,76],
	[47,38,76],
	[9,0,38],
	[23,19,38],

	[127,0,255],
	[191,127,255],
	[82,0,165],
	[124,82,165],
	[63,0,127],
	[95,63,127],
	[38,0,76],
	[57,38,76],
	[19,0,38],
	[28,19,38],

	[191,0,255],
	[223,127,255],
	[124,0,165],
	[145,82,165],
	[95,0,127],
	[111,63,127],
	[57,0,76],
	[66,38,76],
	[28,0,38],
	[33,19,38],

	[255,0,255],
	[255,127,255],
	[165,0,165],
	[165,82,165],
	[127,0,127],
	[127,63,127],
	[76,0,76],
	[76,38,76],
	[38,0,38],
	[38,19,38],

	[255,0,191],
	[255,127,223],
	[165,0,124],
	[165,82,145],
	[127,0,95],
	[127,63,111],
	[76,0,57],
	[76,38,66],
	[38,0,28],
	[38,19,33],

	[255,0,127],
	[255,127,191],
        [165,0,82],
	[165,82,124],
	[127,0,63],
	[127,63,95],
	[76,0,38],
	[76,38,57],
	[38,0,19],
	[38,19,28],

	[255,0,63],
	[255,127,159],
	[165,0,41],
	[165,82,103],
	[127,0,31],
	[127,63,79],
	[76,0,19],
	[76,38,47],
	[38,0,9],
	[38,19,23],

	[0,0,0],
	[45,45,45],
	[91,91,91],
	[137,137,137],
	[183,183,183],
	[179,179,179],
        False                   # bylaer
    ]

    ## function to convert between different color representations 
    def thing2color(self,thing,convert='b',colormap=None,colordict=None):
        def _b2f(c):
            return [ c[0]/255.0,c[1]/255.0,c[2]/255.0 ]
        def _f2b(c):
            return [ round(c[0]*255.0),rount(c[1]*255.0),round(c[2]*255.0) ]
        def _b2i(c):
            for i in range(len(colormap)):
                if colormap[i] == thing:
                    return i
            return False
        def _f2i(c):
            return _b2i(_f2b(c))
        def _isgoodf(x):
            return isinstance(x,float) and x >= 0.0 and x <= 1.0
        def _isgoodb(x):
            return isinstance(x,int) and x >= 0 and x < 256
        def _isgoodfc(c):
            return isinstance(c,(list,tuple)) and len(c) == 3 and\
                _isgoodf(c[0]) and _isgoodf(c[1]) and _isgoodf(c[2])
        def _isgoodbc(c):
            return isinstance(c,(list,tuple)) and len(c) == 3 and\
                _isgoodb(c[0]) and _isgoodb(c[1]) and _isgoodb(c[2])

        if colormap == None:
            colormap = self.colormapAUTOCAD
        if colordict == None:
            colordict = self.colordict
        if convert not in ['f','b','i']:
            raise ValueError('bad colormap conversion')
        c = False
        if _isgoodfc(thing):
            if convert == 'b':
                return _f2b(thing)
            elif convert == 'i':
                return _f2i(thing)
            else:
                return thing
        elif _isgoodbc(thing):
            if convert == 'i':
                return _b2i(thing)
            if convert == 'f':
                return _b2f(thing)
            else:
                return thing
        elif isinstance(thing,str):
            cd = colordict[thing]
            if isinstance(cd,bool) and cd == False:
                raise ValueError('bad color name passed to thing2color: {}'.format(thing))
            if convert == 'i':
                return cd[1]
            else:
                c = cd[0]
        elif isinstance(thing,int):
            if thing < 0 or thing > len(colormap):
                raise ValueError('colormap index out of range: {}'.format(thing))
            if convert == 'i':
                return thing
            else:
                c = colormap[thing]
        else:
            raise ValueError('bad thing passed to thing2color: {}'.format(thing))
        if isinstance(c,bool) and c == False:
            return False
        if convert == 'i':
            return _b2i(c)
        elif convert == 'f':
            return _b2f(c)
        return list(c)
    

