## simple object-oriented framework for dxf-rendered drawable objects
## in yapCAD

import ezdxf
import drawable

## base class to provide dxf drawing functionality
class ezdxfDraw:
    doc = ezdxf.new(dxfversion='R2010')
    doc.layers.new('DEFAULT', dxfattribs={'color': 2})
    msp = doc.modelspace()

    def draw_line(self,p1,p2):
        self.msp.add_line((p1[0], p1[1]), (p2[0], p2[1]))

    def draw_arc(self,p,r,start,end):
        self.msp.add_arc((p[0],p[1]),r,start,end)

    def saveas(self,name):
        self.doc.saveas("{}.dxf".format(name))

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
    
    print("instantiating drawables...")
    print("instantiating Point")
    print("point=Point([10,10])")
    point=Point([10,10],"xo")
    print(point)
    print("instantiating Line")
    print("line=Line([-5,-5],[10,10])")
    line=Line([-5,10],[10,-5])
    print(line)
    print("instantiating Arc")
    print("arc=Arc([0,3],3,45,135)")
    arc=Arc([0,3],3,45,135)
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


    drawable.saveas("test2")
