import pytest
import random
from yapcad.pyglet_drawable import *
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geometry import *
from yapcad.pyglet_drawable import *

"""test functions for the yapcad.geometry module"""


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
        assert L.issampleable()
        assert L.isintersectable()
        assert L.iscontinuous()
        assert not L.isclosed()
        assert close(L.length(),
                     sqrt(10*10+4*4))
        assert vclose(L.center(),
                      center(l))
        
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
        
    def test_visual_intersect(self):
        dd = pygletDraw()
        dd.linecolor='silver'
        # make a polyline spiral
        ply = makeLineSpiral(point(-5,0),2.0,10)
        # wrap it
        Ply = Geometry(ply)
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(5,0),3.0,6.66)
        # wrap it
        Gl = Geometry(gl)
        assert Gl.geom() == gl
        assert Ply.geom() == ply
        # find intersection points
        pts = Gl.intersectXY(Ply)
        # find intersection parameters
        uu = Ply.intersectXY(Gl,params=True)

        dd.draw(Ply)
        dd.draw(Gl)
        dd.polystyle='points'
        dd.pointstyle='x'
        dd.linecolor='aqua'
        dd.draw(pts)
        dd.linecolor='white'
        pts2 = []
        for i in range(len(uu[0])):
            p0 =sample(ply,uu[0][i])
            p1 =sample(gl,uu[1][i])
            assert vclose(p0,p1)
            assert vclose(p0,pts[i])
            pts2.append(p0)
        dd.pointstyle='o'
        dd.draw(pts2)
        dd.display()
