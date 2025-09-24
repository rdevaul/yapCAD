## yapCAD boolen operation support

from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geometry import *
from yapcad.poly import cullZeroLength
class Boolean(Geometry):
    """Boolean operations on Polygons"""

    types = ('union','intersection','difference')

    def __repr__(self):
        return f"Boolean({self.type},{self.elem})"

    def __init__(self,type='union',polys=[]):
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

    @property
    def type(self):
        return self.__type

    def _combine_geom(self,g1,g2):
        gl1 = g1.geom
        gl2 = g2.geom

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
            if len(self.elem)==2:
                self.__outline = self._combine_geom(self.elem[0],self.elem[1])
                self.__outline = cullZeroLength(self.__outline)
                self._setUpdate(False)
                self._setBbox(bbox(self.__outline))
                self._setLength(length(self.__outline))
                self._setCenter(self._calcCenter())
            else:
                raise NotImplementedError(
                    f"don't know how to do {self.type} yet for {len(self.elem)} polygons")
        return deepcopy(self.__outline)
        
        
