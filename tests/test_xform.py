import pytest
from yapcad.xform import *
## unit tests for yapCAD xform.py

class TestXform:
    """unit tests for yapCAD matrix operations"""

    def test_matrix(self):
        foo = Matrix([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
        fooT = Matrix(foo,True)
        bar = Matrix([[1,0,0,1],[0,1,0,1],[0,0,1,1],[0,0,0,1]])
        baz = geom.vect(1,2,3)
        I = Matrix()
        a = 10.0
        assert(I.mul(bar).m == bar.m)
        assert(I.mul(foo).m == foo.m)
        assert(I.mul(fooT).m == fooT.mul(I).m)
        assert(I.mul(I).m == I.m)
        assert(foo.mul(bar).m == [[1,2,3,10],[5,6,7,26],[9,10,11,42],[13,14,15,58]])
        assert(foo.mul(baz) == [18, 46, 74, 102])
        assert(foo.mul(a).m == [[10.0,20.0,30.0,40.0],
                                [50.0,60.0,70.0,80.0],
                                [90.0,100.0,110.0,120.0],
                                [130.0,140.0,150.0,160.0]])
        assert(I.mul(baz) == baz)
        ## homogeneous coordinates test
        assert(geom.homo(foo.mul(baz)) ==
               [18.0/102.0, 46.0/102.0, 74.0/102.0, 1.0])
    
    


