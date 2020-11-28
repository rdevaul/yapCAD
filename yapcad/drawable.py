## base class of drawable for yapCAD

from yapcad.geom import *

## Generic drawing functions -- assumed to use current coordinate
## transform and drawing pen (color, line weight, etc.)
 
class Drawable:
    """simple base class for all drawable geometry"""

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
        if 'x' in self.pointstyle:
            self.draw_x(p,self.pointsize*2)

        if 'o' in self.pointstyle:
            self.draw_circle(p,self.pointsize)

    def __init__(self):
        ## valid pointstyles 'x', 'o', 'xo'
        self.pointstyle = 'xo'
        self.pointsize = 0.1
        self.linestyle = '1'
        self.linecolor = False
        self.fillcolor = False
        ## valid polystyles: 'points', 'lines', 'both'
        self.polystyle = 'lines'
        self.layer = 'default'

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
        # elif isinstance(x,Geometry):
        #     self.draw(x.geom())
        elif isinstance(x,Drawable):
            x.draw()
        else:
            raise ValueError('bad argument to Drawable.draw(): '.format(x))
        
    ## cause drawing page to be rendered -- pure virtual
    def display(self):
        return True

    def set_linecolor(self,c=False):
        if c == False:
            self.linecolor = False
        elif isinstance(c,list) and len(c) == 3 or \
             isinstance(c,int) and c >= 0 or\
             isinstance(c,str) and c in self.colordict:
            self.linecolor = c

    def get_linecolor(self):
        return self.linecolor

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
            if cd == False:
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
        if c == False:
            return False
        if convert == 'i':
            return _b2i(c)
        elif convert == 'f':
            return _b2f(c)
        return list(c)
    

