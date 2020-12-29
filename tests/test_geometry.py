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
        L = Line(a,b)
        assert L.issampleable()
        assert L.isintersectable()
        assert L.iscontinuous()
        assert not L.isclosed()
        assert close(L.length(),
                     sqrt(10*10+4*4))
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
        Ply = Poly(makeLineSpiral(point(-5,0),2.0,10))
        # make an arc-segment geometry list spiral
        Gl = Figure(makeArcSpiral(point(5,0),3.0,6.66,dstep=90.0))
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
            p0 =Gl.sample(uu[1][i])
            assert vclose(p0,pts[i])
            pts2.append(p0)
        dd.pointstyle='o'
        dd.draw(pts2)
        dd.display()

    def test_properties(self):
        # make some randm polys
        box = line(point(-10,-10),
                   point(10,10))
        plys=[]
        nply = 3
        for i in range(nply):
            ply = randomPoly(box,maxr=3.0)
            Ply = Geometry(ply)
            plys.append(Ply)
            assert Ply.issampleable()
            assert Ply.isintersectable()
            assert Ply.iscontinuous()
            assert Ply.isclosed()

        npts = 3
        for i in range(npts):
            ply = randomPoints(box,i+2)
            Ply = Geometry(ply)
            plys.append(Ply)
            assert Ply.issampleable()
            assert Ply.isintersectable()
            assert Ply.iscontinuous()
            assert not Ply.isclosed()

        sprls= 3
        for i in range(sprls):
            r = 3.0
            x = randomCenterInBox(box,r)
            sp = makeArcSpiral(x,r/(i+1.0),i+1.0)
            if not i%2:
                sp.append(line(sample(sp,1.0),
                               sample(sp,0.0)))
            Ply = Geometry(sp)
            plys.append(Ply)
            assert Ply.issampleable()
            assert Ply.isintersectable()
            assert iscontinuousgeomlist(sp)
            assert Ply.iscontinuous()
            if not i%2:
                assert Ply.isclosed()
            else:
                assert not Ply.isclosed()

        dd = pygletDraw()
        dd.linecolor='silver'
        dd.polystyle='both'
        dd.draw(plys)
        dd.display()

    # def test_bool(self):
        
    #     dd = pygletDraw()
    #     dd.linecolor='silver'
    #     dd
        
                        
