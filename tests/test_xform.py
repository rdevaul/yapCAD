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

    def test_translation(self):
        p = geom.vect(1,2,3)
        dx = 4
        dy = -2
        dz = 1
        T = Translation([dx,dy,dz])
        assert(T.mul(p) == geom.vect(1+dx,2+dy,3+dz))
        Tinv = Translation([dx,dy,dz],inverse=True)
        assert(Tinv.mul(p) == geom.vect(1-dx,2-dy,3-dz))
    
    

    def test_rotation_about_a_non_unit_axis(self):
        """Rotation normalises its axis, and the normalising branch used to
        reach for a name that does not exist (geom.geom.scale3), so any axis
        that was not already unit length raised AttributeError. Every axis in
        the suite happened to be unit length, so nothing caught it."""
        import math

        axis = geom.point(1, 1, 0)
        assert not geom.close(geom.mag(axis), 1.0)

        R = Rotation(axis, 45.0)
        unit = geom.point(1 / math.sqrt(2), 1 / math.sqrt(2), 0)
        Runit = Rotation(unit, 45.0)
        for i in range(4):
            for j in range(4):
                assert R.m[i][j] == pytest.approx(Runit.m[i][j], abs=1e-12)

    def test_rotation_rejects_a_zero_length_axis(self):
        with pytest.raises(ValueError):
            Rotation(geom.point(0, 0, 0), 30.0)

