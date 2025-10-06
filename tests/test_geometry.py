import pytest
import os
import random
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geometry import *

"""test functions for the yapcad.geometry module"""

# Control flag for visual tests - set via environment variable or directly
VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw

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
        assert L.sampleable
        assert L.intersectable
        assert L.continuous
        assert not L.closed
        assert close(L.length,
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
        
    @pytest.mark.visual
    def test_visual_intersect(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()
        dd.linecolor='silver'
        # make a polyline spiral
        Ply = Figure(makeLineSpiral(point(-5,0),2.0,10))
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

    @pytest.mark.visual
    def test_properties(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

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

    @pytest.mark.visual
    def test_surface(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()
        dd.linecolor='silver'
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(-10,0),10.0,1.0)
        # "close" the arc with a line segment
        gl.append(line(sample(gl,1.0),
                       sample(gl,0.0)))
        Nautalus = Geometry(gl)
        dd.draw(Nautalus)
        # turn the geometry list into a poly
        surf = Nautalus.surface()
        dd.draw_surface(surf)
        Nautalus.translate(point(0,0,0.5))
        surf2 = Nautalus.surface()
        dd.linecolor='red'
        dd.draw(surf2lines(surf2))
        

        ply2 = [ point(20,20),
                 point(10,20),
                 point(10,-20),
                 point(20,-20),
                 point(20,-10),
                 point(15,-10),
                 point(15,10),
                 point(20,10),
                 point(20,20) ]

        BigC=Geometry(ply2)
        dd.linecolor='silver'
        dd.draw(BigC)
        surf = BigC.surface()
        dd.draw_surface(surf)
        BigC.translate(point(0,0,0.5))
        surf2 = BigC.surface()
        dd.linecolor='red'
        dd.draw(surf2lines(surf2))
        print(f"surface area: {surfacearea(surf)}")
        assert abs(surfacearea(surf)-300.0) < 0.1
        
        # ply2 = translate(ply2,point(0,0,0.5)) # pull it forward
        # dd.linecolor='red'
        # dd.polystyle='both'
        # dd.draw(ply2)
        # surf2,bnd = poly2surface(ply2)
        # dd.draw_surface(surf2)

        # rr1 = makeRoundRect(5,5,0.5,center=point(-15,20))
        # rr2 = makeRoundRect(6,3,0.5,center=point(-12,20))
        # gl3 = Boolean('difference',[rr1,rr2])
        # ply3 = geomlist2poly(gl3.geom)
        # dd.linecolor='aqua'
        # dd.draw(gl3)
        # ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        # surf3,bnd = poly2surface(ply3)
        # dd.linecolor='red'
        # dd.polystyle='both'
        # dd.draw(ply3)
        # dd.draw_surface(surf3)

        # gl4 = [arc(point(0,7),10.0,0.0,270.0),
        #        arc(point(0,7),5.0,0.0,270.0,samplereverse=True)]
        # gl4 = [arc(point(0,7),10.0,0.0,350.0),
        #        arc(point(0,7),5.0,0.0,350.0,samplereverse=True)]
        # dd.linecolor='aqua'
        # dd.draw(gl4)
        # ply4 =geomlist2poly(gl4,minang = 10 )
        # ply4.append(ply4[0])

        # ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        # dd.linecolor='red'
        # dd.polystyle='both'
        # dd.draw(ply4)

        # surf4,bnd = poly2surface(ply4,minlen=1.0)
        # # assert False
        # lines = surf2lines(surf4)
        # lines = translate(lines,point(0,0,0.5))
        # dd.draw_surface(surf4)
        # dd.linecolor='yellow'
        # dd.draw(lines)
        # bnd.append(bnd[0])
        # bnd = translate(bnd,point(0,0,1.0))
        # dd.linecolor='white'
        # dd.draw(bnd)
        
        dd.display()
        
    # def test_bool(self):
        
    #     dd = pygletDraw()
    #     dd.linecolor='silver'
    #     dd
        
                        
