## base class of drawable for yapCAD
## Copyright (c) 2020 Richard W. DeVaul
## Copyright (c) 2020 yapCAD contributors
## All rights reserved
## See licensing terms here: https://github.com/rdevaul/yapCAD/blob/master/LICENSE

from yapcad.geom import *
from yapcad.geometry import *

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

    ## utility function to draw a 2D or 3D bounding box
    def draw_bbox(self,box,dim3=False):
        length =  box[1][0]-box[0][0]
        width = box[1][1]-box[0][1]
        height = box[1][2]-box[0][2]

        p0=box[0]
        p1=add(p0,point(length,0,0))
        p2=add(p1,point(0,width,0))
        p3=add(p0,point(0,width,0))

        self.draw_line(p0,p1)
        self.draw_line(p1,p2)
        self.draw_line(p2,p3)
        self.draw_line(p3,p0)

        if not dim3:
            return
        
        p4=add(p0,point(0,0,height))
        p5=add(p4,point(length,0,0))
        p6=box[1]
        p7=add(p4,point(0,width,0))

        self.draw_line(p4,p5)
        self.draw_line(p5,p6)
        self.draw_line(p6,p7)
        self.draw_line(p7,p4)

        self.draw_line(p0,p4)
        self.draw_line(p1,p5)
        self.draw_line(p2,p6)
        self.draw_line(p3,p7)
        
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
        if pst in ['points','lines','both']:
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
            
        elif iscatmullrom(x):
            from yapcad.spline import sample_catmullrom
            pts = sample_catmullrom(x)
            closed = bool(x[2].get('closed', False))
            for i in range(1, len(pts)):
                self.draw_line(pts[i-1], pts[i])
            if closed and pts:
                self.draw_line(pts[-1], pts[0])
        elif isnurbs(x):
            from yapcad.spline import sample_nurbs
            meta = x[2]
            default = max(32, len(x[1]) * 8)
            count = int(meta.get('samples', default))
            pts = sample_nurbs(x, samples=max(4, count))
            for i in range(1, len(pts)):
                self.draw_line(pts[i-1], pts[i])
            
        elif isinstance(x,Geometry):
            gl = x.geom
            self.draw(gl)
        elif isinstance(x,list): # could be a geometry list, or a list
            # that mixes yapcad.geom elements and yapcad.geometry
            # Geeometry instances.  If the list contains an
            # inappropriate element it will be caught in the "else" case below.
            for e in x:
                self.draw(e)
        else:
            raise ValueError(f'bad argument to Drawable.draw(): {x}')
        
    ## cause drawing page to be rendered -- pure virtual in base class
    def display(self):
        print('pure virtual display function called')
        return True


    ## Standard color names per W3C CSS Color Module Level 4 keywords.
    # Standardized color keywords (CSS Color Module Level 4) plus
    # X11 gray/grey ramp entries for broader palette coverage.
    colordict= {
        'aliceblue': ([240, 248, 255], None),
        'antiquewhite': ([250, 235, 215], None),
        'aqua': ([0, 255, 255], 4),
        'aquamarine': ([127, 255, 212], None),
        'azure': ([240, 255, 255], None),
        'beige': ([245, 245, 220], None),
        'bisque': ([255, 228, 196], None),
        'black': ([0, 0, 0], 0),
        'blanchedalmond': ([255, 235, 205], None),
        'blue': ([0, 0, 255], 5),
        'blueviolet': ([138, 43, 226], None),
        'brown': ([165, 42, 42], None),
        'burlywood': ([222, 184, 135], None),
        'cadetblue': ([95, 158, 160], None),
        'chartreuse': ([127, 255, 0], None),
        'chocolate': ([210, 105, 30], None),
        'coral': ([255, 127, 80], None),
        'cornflowerblue': ([100, 149, 237], None),
        'cornsilk': ([255, 248, 220], None),
        'crimson': ([220, 20, 60], None),
        'cyan': ([0, 255, 255], None),
        'darkblue': ([0, 0, 139], None),
        'darkcyan': ([0, 139, 139], None),
        'darkgoldenrod': ([184, 134, 11], None),
        'darkgray': ([169, 169, 169], None),
        'darkgreen': ([0, 100, 0], None),
        'darkgrey': ([169, 169, 169], None),
        'darkkhaki': ([189, 183, 107], None),
        'darkmagenta': ([139, 0, 139], None),
        'darkolivegreen': ([85, 107, 47], None),
        'darkorange': ([255, 140, 0], None),
        'darkorchid': ([153, 50, 204], None),
        'darkred': ([139, 0, 0], None),
        'darksalmon': ([233, 150, 122], None),
        'darkseagreen': ([143, 188, 143], None),
        'darkslateblue': ([72, 61, 139], None),
        'darkslategray': ([47, 79, 79], None),
        'darkslategrey': ([47, 79, 79], None),
        'darkturquoise': ([0, 206, 209], None),
        'darkviolet': ([148, 0, 211], None),
        'deeppink': ([255, 20, 147], None),
        'deepskyblue': ([0, 191, 255], None),
        'dimgray': ([105, 105, 105], None),
        'dimgrey': ([105, 105, 105], None),
        'dodgerblue': ([30, 144, 255], None),
        'firebrick': ([178, 34, 34], None),
        'floralwhite': ([255, 250, 240], None),
        'forestgreen': ([34, 139, 34], None),
        'fuchsia': ([255, 0, 255], 6),
        'gainsboro': ([220, 220, 220], None),
        'ghostwhite': ([248, 248, 255], None),
        'gold': ([255, 215, 0], None),
        'goldenrod': ([218, 165, 32], None),
        'gray': ([128, 128, 128], 8),
        'gray0': ([0, 0, 0], None),
        'gray1': ([3, 3, 3], None),
        'gray10': ([26, 26, 26], None),
        'gray100': ([255, 255, 255], None),
        'gray11': ([28, 28, 28], None),
        'gray12': ([31, 31, 31], None),
        'gray13': ([33, 33, 33], None),
        'gray14': ([36, 36, 36], None),
        'gray15': ([38, 38, 38], None),
        'gray16': ([41, 41, 41], None),
        'gray17': ([43, 43, 43], None),
        'gray18': ([46, 46, 46], None),
        'gray19': ([48, 48, 48], None),
        'gray2': ([5, 5, 5], None),
        'gray20': ([51, 51, 51], None),
        'gray21': ([54, 54, 54], None),
        'gray22': ([56, 56, 56], None),
        'gray23': ([59, 59, 59], None),
        'gray24': ([61, 61, 61], None),
        'gray25': ([64, 64, 64], None),
        'gray26': ([66, 66, 66], None),
        'gray27': ([69, 69, 69], None),
        'gray28': ([71, 71, 71], None),
        'gray29': ([74, 74, 74], None),
        'gray3': ([8, 8, 8], None),
        'gray30': ([76, 76, 76], None),
        'gray31': ([79, 79, 79], None),
        'gray32': ([82, 82, 82], None),
        'gray33': ([84, 84, 84], None),
        'gray34': ([87, 87, 87], None),
        'gray35': ([89, 89, 89], None),
        'gray36': ([92, 92, 92], None),
        'gray37': ([94, 94, 94], None),
        'gray38': ([97, 97, 97], None),
        'gray39': ([99, 99, 99], None),
        'gray4': ([10, 10, 10], None),
        'gray40': ([102, 102, 102], None),
        'gray41': ([105, 105, 105], None),
        'gray42': ([107, 107, 107], None),
        'gray43': ([110, 110, 110], None),
        'gray44': ([112, 112, 112], None),
        'gray45': ([115, 115, 115], None),
        'gray46': ([117, 117, 117], None),
        'gray47': ([120, 120, 120], None),
        'gray48': ([122, 122, 122], None),
        'gray49': ([125, 125, 125], None),
        'gray5': ([13, 13, 13], None),
        'gray50': ([128, 128, 128], None),
        'gray51': ([130, 130, 130], None),
        'gray52': ([133, 133, 133], None),
        'gray53': ([135, 135, 135], None),
        'gray54': ([138, 138, 138], None),
        'gray55': ([140, 140, 140], None),
        'gray56': ([143, 143, 143], None),
        'gray57': ([145, 145, 145], None),
        'gray58': ([148, 148, 148], None),
        'gray59': ([150, 150, 150], None),
        'gray6': ([15, 15, 15], None),
        'gray60': ([153, 153, 153], None),
        'gray61': ([156, 156, 156], None),
        'gray62': ([158, 158, 158], None),
        'gray63': ([161, 161, 161], None),
        'gray64': ([163, 163, 163], None),
        'gray65': ([166, 166, 166], None),
        'gray66': ([168, 168, 168], None),
        'gray67': ([171, 171, 171], None),
        'gray68': ([173, 173, 173], None),
        'gray69': ([176, 176, 176], None),
        'gray7': ([18, 18, 18], None),
        'gray70': ([178, 178, 178], None),
        'gray71': ([181, 181, 181], None),
        'gray72': ([184, 184, 184], None),
        'gray73': ([186, 186, 186], None),
        'gray74': ([189, 189, 189], None),
        'gray75': ([191, 191, 191], None),
        'gray76': ([194, 194, 194], None),
        'gray77': ([196, 196, 196], None),
        'gray78': ([199, 199, 199], None),
        'gray79': ([201, 201, 201], None),
        'gray8': ([20, 20, 20], None),
        'gray80': ([204, 204, 204], None),
        'gray81': ([207, 207, 207], None),
        'gray82': ([209, 209, 209], None),
        'gray83': ([212, 212, 212], None),
        'gray84': ([214, 214, 214], None),
        'gray85': ([217, 217, 217], None),
        'gray86': ([219, 219, 219], None),
        'gray87': ([222, 222, 222], None),
        'gray88': ([224, 224, 224], None),
        'gray89': ([227, 227, 227], None),
        'gray9': ([23, 23, 23], None),
        'gray90': ([230, 230, 230], None),
        'gray91': ([232, 232, 232], None),
        'gray92': ([235, 235, 235], None),
        'gray93': ([237, 237, 237], None),
        'gray94': ([240, 240, 240], None),
        'gray95': ([242, 242, 242], None),
        'gray96': ([245, 245, 245], None),
        'gray97': ([247, 247, 247], None),
        'gray98': ([250, 250, 250], None),
        'gray99': ([252, 252, 252], None),
        'green': ([0, 127, 0], 94),
        'greenyellow': ([173, 255, 47], None),
        'grey': ([128, 128, 128], None),
        'grey0': ([0, 0, 0], None),
        'grey1': ([3, 3, 3], None),
        'grey10': ([26, 26, 26], None),
        'grey100': ([255, 255, 255], None),
        'grey11': ([28, 28, 28], None),
        'grey12': ([31, 31, 31], None),
        'grey13': ([33, 33, 33], None),
        'grey14': ([36, 36, 36], None),
        'grey15': ([38, 38, 38], None),
        'grey16': ([41, 41, 41], None),
        'grey17': ([43, 43, 43], None),
        'grey18': ([46, 46, 46], None),
        'grey19': ([48, 48, 48], None),
        'grey2': ([5, 5, 5], None),
        'grey20': ([51, 51, 51], None),
        'grey21': ([54, 54, 54], None),
        'grey22': ([56, 56, 56], None),
        'grey23': ([59, 59, 59], None),
        'grey24': ([61, 61, 61], None),
        'grey25': ([64, 64, 64], None),
        'grey26': ([66, 66, 66], None),
        'grey27': ([69, 69, 69], None),
        'grey28': ([71, 71, 71], None),
        'grey29': ([74, 74, 74], None),
        'grey3': ([8, 8, 8], None),
        'grey30': ([76, 76, 76], None),
        'grey31': ([79, 79, 79], None),
        'grey32': ([82, 82, 82], None),
        'grey33': ([84, 84, 84], None),
        'grey34': ([87, 87, 87], None),
        'grey35': ([89, 89, 89], None),
        'grey36': ([92, 92, 92], None),
        'grey37': ([94, 94, 94], None),
        'grey38': ([97, 97, 97], None),
        'grey39': ([99, 99, 99], None),
        'grey4': ([10, 10, 10], None),
        'grey40': ([102, 102, 102], None),
        'grey41': ([105, 105, 105], None),
        'grey42': ([107, 107, 107], None),
        'grey43': ([110, 110, 110], None),
        'grey44': ([112, 112, 112], None),
        'grey45': ([115, 115, 115], None),
        'grey46': ([117, 117, 117], None),
        'grey47': ([120, 120, 120], None),
        'grey48': ([122, 122, 122], None),
        'grey49': ([125, 125, 125], None),
        'grey5': ([13, 13, 13], None),
        'grey50': ([128, 128, 128], None),
        'grey51': ([130, 130, 130], None),
        'grey52': ([133, 133, 133], None),
        'grey53': ([135, 135, 135], None),
        'grey54': ([138, 138, 138], None),
        'grey55': ([140, 140, 140], None),
        'grey56': ([143, 143, 143], None),
        'grey57': ([145, 145, 145], None),
        'grey58': ([148, 148, 148], None),
        'grey59': ([150, 150, 150], None),
        'grey6': ([15, 15, 15], None),
        'grey60': ([153, 153, 153], None),
        'grey61': ([156, 156, 156], None),
        'grey62': ([158, 158, 158], None),
        'grey63': ([161, 161, 161], None),
        'grey64': ([163, 163, 163], None),
        'grey65': ([166, 166, 166], None),
        'grey66': ([168, 168, 168], None),
        'grey67': ([171, 171, 171], None),
        'grey68': ([173, 173, 173], None),
        'grey69': ([176, 176, 176], None),
        'grey7': ([18, 18, 18], None),
        'grey70': ([178, 178, 178], None),
        'grey71': ([181, 181, 181], None),
        'grey72': ([184, 184, 184], None),
        'grey73': ([186, 186, 186], None),
        'grey74': ([189, 189, 189], None),
        'grey75': ([191, 191, 191], None),
        'grey76': ([194, 194, 194], None),
        'grey77': ([196, 196, 196], None),
        'grey78': ([199, 199, 199], None),
        'grey79': ([201, 201, 201], None),
        'grey8': ([20, 20, 20], None),
        'grey80': ([204, 204, 204], None),
        'grey81': ([207, 207, 207], None),
        'grey82': ([209, 209, 209], None),
        'grey83': ([212, 212, 212], None),
        'grey84': ([214, 214, 214], None),
        'grey85': ([217, 217, 217], None),
        'grey86': ([219, 219, 219], None),
        'grey87': ([222, 222, 222], None),
        'grey88': ([224, 224, 224], None),
        'grey89': ([227, 227, 227], None),
        'grey9': ([23, 23, 23], None),
        'grey90': ([230, 230, 230], None),
        'grey91': ([232, 232, 232], None),
        'grey92': ([235, 235, 235], None),
        'grey93': ([237, 237, 237], None),
        'grey94': ([240, 240, 240], None),
        'grey95': ([242, 242, 242], None),
        'grey96': ([245, 245, 245], None),
        'grey97': ([247, 247, 247], None),
        'grey98': ([250, 250, 250], None),
        'grey99': ([252, 252, 252], None),
        'honeydew': ([240, 255, 240], None),
        'hotpink': ([255, 105, 180], None),
        'indianred': ([205, 92, 92], None),
        'indigo': ([75, 0, 130], None),
        'ivory': ([255, 255, 240], None),
        'khaki': ([240, 230, 140], None),
        'lavender': ([230, 230, 250], None),
        'lavenderblush': ([255, 240, 245], None),
        'lawngreen': ([124, 252, 0], None),
        'lemonchiffon': ([255, 250, 205], None),
        'lightblue': ([173, 216, 230], None),
        'lightcoral': ([240, 128, 128], None),
        'lightcyan': ([224, 255, 255], None),
        'lightgoldenrodyellow': ([250, 250, 210], None),
        'lightgray': ([211, 211, 211], None),
        'lightgreen': ([144, 238, 144], None),
        'lightgrey': ([211, 211, 211], None),
        'lightpink': ([255, 182, 193], None),
        'lightsalmon': ([255, 160, 122], None),
        'lightseagreen': ([32, 178, 170], None),
        'lightskyblue': ([135, 206, 250], None),
        'lightslategray': ([119, 136, 153], None),
        'lightslategrey': ([119, 136, 153], None),
        'lightsteelblue': ([176, 196, 222], None),
        'lightyellow': ([255, 255, 224], None),
        'lime': ([0, 255, 0], 3),
        'limegreen': ([50, 205, 50], None),
        'linen': ([250, 240, 230], None),
        'magenta': ([255, 0, 255], None),
        'maroon': ([127, 0, 0], 14),
        'mediumaquamarine': ([102, 205, 170], None),
        'mediumblue': ([0, 0, 205], None),
        'mediumorchid': ([186, 85, 211], None),
        'mediumpurple': ([147, 112, 219], None),
        'mediumseagreen': ([60, 179, 113], None),
        'mediumslateblue': ([123, 104, 238], None),
        'mediumspringgreen': ([0, 250, 154], None),
        'mediumturquoise': ([72, 209, 204], None),
        'mediumvioletred': ([199, 21, 133], None),
        'midnightblue': ([25, 25, 112], None),
        'mintcream': ([245, 255, 250], None),
        'mistyrose': ([255, 228, 225], None),
        'moccasin': ([255, 228, 181], None),
        'navajowhite': ([255, 222, 173], None),
        'navy': ([0, 0, 128], 174),
        'oldlace': ([253, 245, 230], None),
        'olive': ([127, 127, 0], 54),
        'olivedrab': ([107, 142, 35], None),
        'orange': ([255, 165, 0], None),
        'orangered': ([255, 69, 0], None),
        'orchid': ([218, 112, 214], None),
        'palegoldenrod': ([238, 232, 170], None),
        'palegreen': ([152, 251, 152], None),
        'paleturquoise': ([175, 238, 238], None),
        'palevioletred': ([219, 112, 147], None),
        'papayawhip': ([255, 239, 213], None),
        'peachpuff': ([255, 218, 185], None),
        'peru': ([205, 133, 63], None),
        'pink': ([255, 192, 203], None),
        'plum': ([221, 160, 221], None),
        'powderblue': ([176, 224, 230], None),
        'purple': ([128, 0, 128], 214),
        'rebeccapurple': ([102, 51, 153], None),
        'red': ([255, 0, 0], 1),
        'rosybrown': ([188, 143, 143], None),
        'royalblue': ([65, 105, 225], None),
        'saddlebrown': ([139, 69, 19], None),
        'salmon': ([250, 128, 114], None),
        'sandybrown': ([244, 164, 96], None),
        'seagreen': ([46, 139, 87], None),
        'seashell': ([255, 245, 238], None),
        'sienna': ([160, 82, 45], None),
        'silver': ([192, 192, 192], 9),
        'skyblue': ([135, 206, 235], None),
        'slateblue': ([106, 90, 205], None),
        'slategray': ([112, 128, 144], None),
        'slategrey': ([112, 128, 144], None),
        'snow': ([255, 250, 250], None),
        'springgreen': ([0, 255, 127], None),
        'steelblue': ([70, 130, 180], None),
        'tan': ([210, 180, 140], None),
        'teal': ([0, 127, 127], 134),
        'thistle': ([216, 191, 216], None),
        'tomato': ([255, 99, 71], None),
        'turquoise': ([64, 224, 208], None),
        'violet': ([238, 130, 238], None),
        'wheat': ([245, 222, 179], None),
        'white': ([255, 255, 255], 7),
        'whitesmoke': ([245, 245, 245], None),
        'yellow': ([255, 255, 0], 2),
        'yellowgreen': ([154, 205, 50], None),
    }

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
            return [c[0] / 255.0, c[1] / 255.0, c[2] / 255.0]

        def _f2b(c):
            return [round(c[0] * 255.0), round(c[1] * 255.0), round(c[2] * 255.0)]

        def _nearest_index(rgb):
            best_idx = 0
            best_dist = float('inf')
            for idx, candidate in enumerate(colormap):
                if not isinstance(candidate, (list, tuple)) or len(candidate) != 3:
                    continue
                dist = ((candidate[0] - rgb[0]) ** 2 +
                        (candidate[1] - rgb[1]) ** 2 +
                        (candidate[2] - rgb[2]) ** 2)
                if dist < best_dist:
                    best_idx = idx
                    best_dist = dist
                    if dist == 0:
                        break
            return best_idx

        def _b2i(c):
            for i in range(len(colormap)):
                if colormap[i] == c:
                    return i
            return _nearest_index(c)

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
            key = thing.lower()
            if key not in colordict:
                raise ValueError('bad color name passed to thing2color: {}'.format(thing))
            cd = colordict[key]
            if isinstance(cd,bool) and cd == False:
                raise ValueError('bad color name passed to thing2color: {}'.format(thing))
            if convert == 'i':
                idx = cd[1]
                if idx is None:
                    idx = _nearest_index(cd[0])
                return idx
            else:
                c = list(cd[0])
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
    
