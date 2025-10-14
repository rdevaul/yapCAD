## yapCAD boolen operation support for 2D closed curves. For three dimensional
## boolean operations, see geom3d.py

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geometry import *
from yapcad.poly import cullZeroLength
class Boolean(Geometry):
    """Boolean operations on Polygons"""

    types = ('union','intersection','difference')

    def __repr__(self):
        return f"Boolean({self.type},{self.elem})"

    def __init__(self,type='union',polys=[], *, minang=5.0, minlen=0.25):
        super().__init__()
        self._setClosed(True)
        self._setSampleable(True)
        if not type in self.types:
            raise ValueError('invalid type passed to Boolean(): {}'.format(type))
        for p in polys:
            if not ( isinstance(p,Geometry) and p.isclosed()):
                raise ValueError('not closed Geometry instance: {}'.format(p))
            self.elem.append(deepcopy(p))

        self.__type=type
        self._setUpdate(True)
        self.__outline=[]
        self.__minang = float(minang)
        self.__minlen = float(minlen)

    @property
    def type(self):
        return self.__type

    def _poly_to_segments(self, poly):
        segments = []
        if not poly:
            return segments
        for i in range(1, len(poly)):
            segments.append([poly[i - 1], poly[i]])
        return segments

    def _prepare_geom(self, geom_obj):
        if isinstance(geom_obj, Geometry):
            gl = geom_obj.geom
        else:
            gl = geom_obj

        # If gl is a single arc/line/poly, wrap it in a list for geomlist2poly_with_holes
        if isarc(gl) or isline(gl) or ispoly(gl):
            gl = [gl]

        try:
            outer, holes = geomlist2poly_with_holes(gl, self.__minang, self.__minlen, checkcont=False)
        except Exception:
            if isgeomlist(gl):
                return gl
            return []

        if not outer:
            return []

        segments = self._poly_to_segments(outer)
        for hole in holes:
            segments.extend(self._poly_to_segments(hole))
        return segments

    def _combine_geom(self,g1,g2):
        gl1 = self._prepare_geom(g1)
        gl2 = self._prepare_geom(g2)

        return combineglist(gl1,gl2,self.type)

    def _calcCenter(self):
        l = len(self.__outline)
        if l == 0:
            return None
        elif l == 1:
            return center(self.__outline[0]) # center is sole point
        else:
            
            if dist(center(self.__outline[0]),
                    center(self.__outline[-1])) < epsilon:
                l -= 1
            
            p = center(self.__outline[0])
            for i in range(1,l):
                p = add(center(self.__outline[i]),p)

            return scale3(p,1/l)

    def translate(self,delta):
        self._setUpdate(True)
        for p in self.elem:
            p.translate(delta)

    def scale(self,sx,sy=False,sz=False,cent=point(0,0)):
        self._setUpdate(True)
        for p in self.elem:
            p.scale(sx,sy,sz,cent)

    def rotate(self,angle,cent=point(0,0,0),axis=point(0,0,1)):
        self._setUpdate(True)
        for p in self.elem:
            p.rotate(angle,cent,axis)

    def mirror(self,plane):
        self._setUpdate(True)
        for p in self.elem:
            p.mirror(plane)

    def transform(self,m):
        self._setUpdate(True)
        for p in self.elem:
            p.transform(m)
            
    @property
    def geom(self):
        if self.update:
            if len(self.elem) < 2:
                raise ValueError('Boolean requires at least two operands')

            result = self.elem[0]
            for operand in self.elem[1:]:
                result = self._combine_geom(result, operand)

            self.__outline = cullZeroLength(result)
            self._setUpdate(False)
            self._setBbox(bbox(self.__outline))
            self._setLength(length(self.__outline))
            self._setCenter(self._calcCenter())
        return deepcopy(self.__outline)
        
        
