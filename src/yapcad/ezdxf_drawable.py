## simple yapCAD framework for dxf-rendered drawable objects using
## ezdxf package.
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

from yapcad.geom import *
import yapcad.drawable as drawable
import ezdxf
from ezdxf.enums import TextEntityAlignment

## class to provide dxf drawing functionality
class ezdxfDraw(drawable.Drawable):

    def __init__(self):
        super().__init__()
        
        # setup=False avoids creating default blocks (like _CLOSEDFILLED) that
        # contain SOLID entities unsupported by some CAD programs (e.g., FreeCAD)
        self.__doc = ezdxf.new(dxfversion='R2010', setup=False)
        self.__doc.header['$MEASUREMENT'] = 1 # metric
        self.__doc.header['$INSUNITS'] = 4 # millimeters
        self.__doc.layers.new('PATHS',  dxfattribs={'color': 7}) #white
        self.__doc.layers.new('DRILLS',  dxfattribs={'color': 4}) #aqua
        self.__doc.layers.new('DOCUMENTATION', dxfattribs={'color': 2}) #yellow
        self.__linetypelist = [ False,'Continuous']
        self.__msp = self.__doc.modelspace()
        self.__filename = "yapCAD-out"
        self.layerlist = [False, '0','PATHS','DRILLS','DOCUMENTATION']


    def __repr__(self):
        return 'an instance of ezdxfDraw'

    ## properties
    
    # @drawable.Drawable.layer.setter
    # def layer(self,layer=False):
    #     if not layer in self.layerlist:
    #         raise ValueError('bad layer passed to layerset: {}'.format(layer))
    #     self._set_layer(layer)

    @drawable.Drawable.linetype.setter
    def linetype(self,linetype=False):
        if not linetype in self.__linetypelist:
            raise ValueError('bad linetype in linetypeset: {}'.format(linetype))
        self._set_linetype(linetype)

    @property
    def filename(self):
        return self.__filename

    def _set_filename(self,name):
        self.__filename = name

    @filename.setter
    def filename(self,name):
        if not isinstance(name,str):
            raise ValueError('bad (non-string) filename: '+str(name))
        self._set_filename(name)
        
    ## Overload virtual yapcasd.drawable base class drawing methods
    
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
            
        self.__msp.add_line((p1[0], p1[1]), (p2[0], p2[1]),
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
            self.__msp.add_circle((p[0],p[1]),r,
                          dxfattribs={'layer': layer,
                                      'color': color,
                                      'linetype': linetype})
        else:
            self.__msp.add_arc((p[0],p[1]),r,start,end,
                          dxfattribs={'layer': layer,
                                      'color': color,
                                      'linetype': linetype})

    def draw_ellipse(self, center, semi_major, semi_minor, rotation, start, end):
        """Draw an ellipse or elliptical arc to DXF.

        DXF ellipse is defined by:
        - center point
        - major axis endpoint relative to center
        - ratio of minor to major axis
        - start/end parameters (0 to 2*pi)
        """
        import math

        layer = self.layer
        if layer == False:
            layer = '0'
        color = self.linecolor
        if color:
            color = self.thing2color(self.linecolor, 'i')
        else:
            color = 256  # bylayer
        linetype = self.linetype
        if linetype == False:
            linetype = 'Continuous'

        # Convert rotation from degrees to radians
        rot_rad = math.radians(rotation)

        # Major axis endpoint relative to center
        major_axis = (semi_major * math.cos(rot_rad),
                      semi_major * math.sin(rot_rad),
                      0)

        # Ratio of minor to major
        ratio = semi_minor / semi_major

        # Convert start/end angles from degrees to radians
        # DXF ellipse uses parametric angles (0 to 2*pi)
        if start == 0 and end == 360:
            start_param = 0.0
            end_param = math.tau
        else:
            start_param = math.radians(start)
            end_param = math.radians(end)
            if end_param < start_param:
                end_param += math.tau

        self.__msp.add_ellipse(
            center=(center[0], center[1]),
            major_axis=major_axis,
            ratio=ratio,
            start_param=start_param,
            end_param=end_param,
            dxfattribs={'layer': layer,
                        'color': color,
                        'linetype': linetype}
        )

    def draw_text(self,text,location,
                  align=TextEntityAlignment.LEFT,
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

        alignment = [];
        if align == 'LEFT':
            alignment = TextEntityAlignment.LEFT
        elif align == 'CENTER':
            alignment = TextEntityAlignment.CENTER
        elif align == 'RIGHT':
            alignment = TextEntityAlignment.RIGHT
        else:
            ## default alignment is left
            alignment = TextEntityAlignment.LEFT
        
        dxfattr = dict()
        if 'style' in attr:
            dxfattr['style'] = attr['style']
        if 'height' in attr:
            dxfattr['height'] = attr['height']
        dxfattr.update({ 'layer': layer,
                         'color': color})
        self.__msp.add_text(
            text,
            dxfattribs=dxfattr).set_placement(
                (location[0],location[1]),
                align=alignment
            )

    def display(self):
        self.__doc.saveas("{}.dxf".format(self.filename))

