import pytest
from geom import *
## unit tests for yapCAD geom.py

class TestPoint:
    """unit tests for yapCAD point functions"""

    def test_create(self):
        a = point(5,0)
        b = point(0,5,-2)
        c = point(-2.3,4.6,-9.2,0.5)
        bb = point(b)
        assert a == [5,0,0,1]
        assert b == [0,5,-2,1]
        assert c == [-2.3,4.6,-9.2,0.5]
        assert bb[0] == b[0] and bb[1] == b[1] \
            and bb[2] == b[2] and bb[3] == b[3]

    def test_discrimate(self):
        a = point(5,0)
        b = point(0,5)
        assert ispoint(a)
        assert ispoint(b)
        assert ispoint([0,2,2,1])
        assert not ispoint([1,2])

    def test_format(self):
        a = point(5,0)
        b = point(2,3,2)
        c = point(1,2,3,4)
        assert vstr(a) == '[5, 0]'
        assert vstr(b) == '[2, 3, 2]'
        assert vstr(c) == '[1, 2, 3, 4]'

class TestLine:
    def test_sample(self):
        a = point(-5,-1)
        b = point(5,3)
        l = line(a,b)
        p0 =sampleline(l,-0.5)
        p1 =sampleline(l,0.0)
        p2 =sampleline(l,0.5)
        p3 =sampleline(l,1.0)
        p4 =sampleline(l,1.5)
        assert mag(sub(p0,point(-10,-3))) < epsilon
        assert mag(sub(p1,a)) < epsilon
        assert mag(sub(p2,point(0,1))) < epsilon
        assert mag(sub(p3,b)) < epsilon
        assert mag(sub(p4,point(10,5))) < epsilon
        u0 = unsampleline(l,p0)
        u1 = unsampleline(l,p1)
        u2 = unsampleline(l,p2)
        u3 = unsampleline(l,p3)
        u4 = unsampleline(l,p4)
        assert abs(u0+0.5) < epsilon
        assert abs(u1) < epsilon
        assert abs(u2-0.5) < epsilon
        assert abs(u3-1.0) < epsilon
        assert abs(u4-1.5) < epsilon
        p = point(100,100)
        assert not unsampleline(l,p)
        
class TestPoly:
    def test_create(self):
        with pytest.raises(ValueError):
            poly()

        a = point(0,-5)
        b = point(5,0)
        c = point(0,5)
        d = point(-5,0)

        pol1 = poly(a,b,c)
        assert pol1 == [a,b,c]

        pol2 = poly(a,b,c,d)
        assert pol2 == [a,b,c,d]
        
        pol3 = poly(pol1)
        assert not pol2 == pol1
        assert pol3 == pol1

        with pytest.raises(ValueError):
            poly(a,b,'foo')

        with pytest.raises(ValueError):
            poly(a,b)

    def test_discriminate(self):
        a = point(0,-5)
        b = point(5,0)
        c = point(0,5)
        
        pol1 = poly(a,b,c)
        pol2 = poly(a,b,c,a)
        
        assert ispoly(pol1)
        assert not ispoly(a)
        assert not ispoly('foo')

        assert ispoly(pol2)
        assert ispolygon(pol2)
        assert not ispolygon(pol1)
        assert not ispolygon([])
