## yapCAD geometry-generating superclass
## =====================================

## subclasses of Geometry implement the geom() method, which produce
## geometry lists or computational geometry primitives, as described
## by the functions in geom.py

## subclasses of SampleGeometry implelement geom() and sample()
## methods.  The geom() method is guaranteed to produce a geometry
## list, and the sample() method allows for parametric sampling, which
## is guaranteed to "make sense" in the interval 0.0 <= u <= 1.0.
## Unless there is a good reason otherwise, subclasses of
## SampleGeometry should guarantee C0 continuity*

## * For fractal elements, or other complex forms where C0 continuity
## is less meaningful, SampleGeometry might still be approprite.
## However, if you have diconnected line and arc segments, with a
## classically-definalble gap (or the possibility thereof) then
## implement Geometry instead

import copy
from yapcad.geom import *


class Geometry:
    """ generalized computational geometry class """

    def __repr__(self):
        return 'geometry base class wrapper for: {}'.format(vstr(self._elem))

    def __init__(self,a=False):
        self._update=True
        self._elem=[]
        if a:
            if ispoint(a) or isline (a) or isarc(a) or ispoly(a) \
                or isgeomlist(a) or isinstance(a,Geometry):
                self._elem=[ deepcopy(a) ]
            else:
                raise ValueError('bad argument to Geometry class constructor: {}'.format(a))
        

    def _updateInternals(self):
        return

    def geom(self):
        if self.update:
            self._updateInternals()
        return deepcopy(self._elem)

    

    
class SampleGeometry(Geometry):
    """ generalized sampleable geometry class"""

    def __init__(self,a=False):
        super().__init__(a)

    def __repr__(self):
        return 'sampleable geometry base class wrapper for: {}'.format(vstr(self._elem))
        
    def sample(self,u):
        if self._update:
            self._updateInternals()
        return sample(self.geom(),u)
    

class IntersectGeometry(SampleGeometry):
    """ generalized intersectable geometry class"""

    def __init__(self,a=False):
        super().__init__(a)

    def __repr__(self):
        return 'intersectable geometry base class wrapper for: {}'.format(vstr(self._elem))
        
    def intersectXY(self,g,inside=True,params=False):
        if self._update:
            self._updateInternals()
        return intersectXY(g,self.geom(),inside,params)
    
