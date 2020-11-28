## simple yapCAD framework for dxf-rendered drawable objects using
## ezdxf package.
## Original author: Richard W. DeVaul

from yapcad.geom import *
import yapcad.drawable as drawable
import ezdxf

## class to provide dxf drawing functionality
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
        self.linecolor=False    # default is 256 = autocad BYLAYER
        self.linetype=False # default is 'Continuous'
        
        super().__init__()

    def __repr__(self):
        return 'an instance of ezdxfDraw'

    def set_layer(self,layer=False):
        if not layer in self.layerlist:
            raise ValueError('bad layer passed to layerset: {}'.format(layer))
        self.layer=layer
        
    def set_linetype(self,linetype=False):
        if not linetype in self.linetypelist:
            raise ValueError('bad linetype in linetypeset: {}'.format(linetype))
        self.linetype=linetype
        
    def draw_line(self,p1,p2):
        layer=self.layer
        if layer == False:
            layer = '0'
        color = self.linecolor
        if color:
            color = self.thing2color(self.linecolor,'i')
        else:
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
        color = self.linecolor
        if color:
            color = self.thing2color(self.linecolor,'i')
        else:
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

        color = []
        if 'color' in attr:
            color = self.thing2color(attr['color'],'i')
        elif self.linecolor:
            color = self.thing2color(self.linecolor,'i')
        else:
            color = 256 # by layer

        dxfattr = dict()
        if 'style' in attr:
            dxfattr['style'] = attr['style']
        if 'height' in attr:
            dxfattr['height'] = attr['height']
        dxfattr.update({ 'layer': layer,
                         'color': color})
        self.msp.add_text(text,dxfattr).set_pos((location[0],location[1]),
                                             align=align)
            
    def set_filename(self,name):
        self.filename = name
        
    def display(self):
        self.doc.saveas("{}.dxf".format(self.filename))

