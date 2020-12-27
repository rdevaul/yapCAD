import pytest
import random
from yapcad.geom import *
from yapcad.geometry import *
from yapcad.pyglet_drawable import *

"""test functions for the yapcad.geometry module"""

def makeSpiral(center, turnRad, # radius after one full turn
               turns, # number of turns
               dstep = 5.0): # sampling resolution in degrees
    # make a spiral of points
    spiral = []
    rstep = turnRad*dstep/360.0

    for i in range(round(360*turns/3)):
        ang = i * dstep*pi2/360.0
        r = i * rstep
        p = add(center,
                point(math.cos(ang)*r,math.sin(ang)*r))
        spiral.append(p)
    return spiral


class TestGeometry:
    """Tests for the Geometry class"""
    ## make some test geometry

    # bounding box for tests
    bb1 = [point(-10,-10,0),
           point(10,10,0)]

    def test_geom(self):
        """test basic yapcad.geom figure wrapping"""

        a = point(-5,-1)
        b = point(5,3)
        l = line(a,b)
        L = Geometry(l)
        
        print("line sample tests")
        p0 =L.sample(-0.5)
        p1 =L.sample(0.0)
        p2 =L.sample(0.5)
        p3 =L.sample(1.0)
        p4 =L.sample(1.5)
        assert vclose(p0,point(-10,-3))
        assert vclose(p1,a)
        assert vclose(p2,point(0,1))
        assert vclose(p3,b)
        assert vclose(p4,point(10,5))
        print("line unsample tests")
        u0 = L.unsample(p0)
        u1 = L.unsample(p1)
        u2 = L.unsample(p2)
        u3 = L.unsample(p3)
        u4 = L.unsample(p4)
        assert close(u0,-0.5)
        assert close(u1,0.0)
        assert close(u2,0.5)
        assert close(u3,1.0)
        assert close(u4,1.5)
        p = point(100,100)
        assert not L.unsample(p)
    
    
    
