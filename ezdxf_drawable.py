## simple object-oriented framework for dxf-rendered drawable objects
## in yapCAD

from geom import *
import ezdxf
import drawable

## base class to provide dxf drawing functionality
class ezdxfDraw(drawable.Drawable):

    def __init__(self):
        self.doc = ezdxf.new(dxfversion='R2010',setup=True)
        self.doc.layers.new('PATHS',  dxfattribs={'color': 7}) #white
        self.doc.layers.new('DRILLS',  dxfattribs={'color': 4}) #aqua
        self.doc.layers.new('DOCUMENTATION', dxfattribs={'color': 2}) #yellow
        self.layerlist = [False, '0','PATHS','DRILLS','DOCUMENTATION']
        self.linetypelist = [ False,'Continuous']
        self.msp = self.doc.modelspace()
        self.filename = "yapCAD-out"

        self.layer=False        # no layer set, use '0'
        self.color=False    # default is 256 = autocad BYLAYER
        self.linetype=False # default is 'Continuous'
        
        super().__init__()

    def __repr__(self):
        return 'an instance of ezdxfDraw'

    def layerset(self,layer=False):
        if not layer in self.layerlist:
            raise ValueError('bad layer passed to layerset: {}'.format(layer))
        self.layer=layer
        
    def colorset(self,color=False):
        if not color in range(257):
            raise ValueError('bad color in colorset: {}'.format(color))
        self.color=color

    def linetypeset(self,linetype=False):
        if not linetype in self.linetypelist:
            raise ValueError('bad linetype in linetypeset: {}'.format(linetype))
        self.linetype=linetype
        
    def draw_line(self,p1,p2):
        layer=self.layer
        if layer == False:
            layer = '0'
        color = self.color
        if color == False:
            color = 256 # bylaer
        linetype = self.linetype
        if linetype == False:
            linetype = 'Continuous'
            
        self.msp.add_line((p1[0], p1[1]), (p2[0], p2[1]),
                          dxfattribs={'layer': layer,
                                      'color': color,
                                      'linetype': linetype})
            

    def draw_arc(self,p,r,start,end):
        layer=self.layer
        if layer == False:
            layer = '0'
        color = self.color
        if color == False:
            color = 256 # bylaer
        linetype = self.linetype
        if linetype == False:
            linetype = 'Continuous'

        if start==0 and end==360:
            self.msp.add_circle((p[0],p[1]),r,
                          dxfattribs={'layer': layer,
                                      'color': color,
                                      'linetype': linetype})
        else:
            self.msp.add_arc((p[0],p[1]),r,start,end,
                          dxfattribs={'layer': layer,
                                      'color': color,
                                      'linetype': linetype})

    def draw_text(self,text,location,
                  align='LEFT',
                  attr={'style': 'LiberationMono',
                        'height': .75}):
        layer=self.layer
        if layer == False:
            layer = '0'
        color = self.color
        if color == False:
            color = 256 # bylaer
        dxfattr = dict(attr)
        dxfattr.update({ 'layer': layer,
                         'color': color})
        self.msp.add_text(text,dxfattr).set_pos((location[0],location[1]),
                                             align=align)
            
    def saveas(self,name):
        self.filename = name
        
    def display(self):
        self.doc.saveas("{}.dxf".format(self.filename))

## multiply-inherited drawing classes
class Point(ezdxfDraw,drawable.Point):
    """ezdxf-rendering Point class"""
        
class Line(ezdxfDraw,drawable.Line):
    """ezdxf-rendering Line class"""

class Arc(ezdxfDraw,drawable.Arc):
    """ezdxf-rendering Arc class"""

if __name__ == "__main__":
    print("testing for drawable.py")
    print("-----------------------")
    print("instantiating ezdxfDraw")
    drawable=ezdxfDraw()
    print("setting the save-as filename")
    drawable.saveas("test2")
    
    print("instantiating drawables...")
    print("instantiating Point")
    print("point=Point(vect(10,10))")
    point=Point(vect(10,10),"xo")
    print(point)
    print("instantiating Line")
    print("line=Line(vect(-5,-5),vect(10,10))")
    line=Line(vect(-5,10),vect(10,-5))
    print(line)
    print("instantiating Arc")
    print("arc=Arc(vect(0,3),3,45,135)")
    arc=Arc(vect(0,3),3,45,135)
    print(arc)
    print("draw tests")
    point.draw()
    line.draw()
    arc.draw()

    print("sample test")
    print("drawing 4 samples from line")
    line.draw()
    for i in range(4):
        u = i/3
        p = line.sample(u)
        pnt = Point(p,"o")
        pnt.draw()
    print("drawing 8 samples from arc")
    arc.draw()
    for i in range(8):
        u = i/7
        p = arc.sample(u)
        pnt = Point(p,"o")
        pnt.draw()

    ## render results
    drawable.display()
