import pytest
from yapcad.geom import *
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
        assert not ispoint(vect(1,2,3,-1))
        assert ispoint([0,2,2,1])
        assert not ispoint([1,2])

    def test_format(self):
        a = point(5,0)
        b = point(2,3,2)
        c = point(1,2,3,4)
        assert vstr(a) == '[5, 0]'
        assert vstr(b) == '[2, 3, 2]'
        assert vstr(c) == '[1, 2, 3, 4]'

class TestOperations:
    def test_vect(self):
        
        a = point(5,0)
        b = point(0,5)
        c = point(-3,-3)
        d = point(1,1)
        e = point(10,10)
        assert close(mag(a),5.0)
        assert vclose(add(a,b),point(5,5))
        assert vclose(sub(a,b),point(5,-5))
        assert close(dot(a,b),0)
        assert close(dot(d,c),-6)
        assert vclose(cross(a,b),point(0,0,25))
        assert vclose(cross(b,a),point(0,0,-25))
        assert close(mag(sub(a,b)),sqrt(50))

class TestLine:
    def test_create(self):
        
        a = point(5,0)
        b = point(0,5)
        c = point(-3,-3)
        d = point(1,1)
        e = point(10,10)
        f = vect(1,2,3,-1)

        assert line(a,b) == [a,b]
        assert isline(line(a,b))
        assert isline([a,b])
        assert not isline([a,f])
        with pytest.raises(ValueError):
            assert not lsine(line(a,f))
        
    def test_sample(self):
        a = point(-5,-1)
        b = point(5,3)
        l = line(a,b)
        print("sample tests")
        p0 =sampleline(l,-0.5)
        p1 =sampleline(l,0.0)
        p2 =sampleline(l,0.5)
        p3 =sampleline(l,1.0)
        p4 =sampleline(l,1.5)
        assert vclose(p0,point(-10,-3))
        assert vclose(p1,a)
        assert vclose(p2,point(0,1))
        assert vclose(p3,b)
        assert vclose(p4,point(10,5))
        print("unsample tests")
        u0 = unsampleline(l,p0)
        u1 = unsampleline(l,p1)
        u2 = unsampleline(l,p2)
        u3 = unsampleline(l,p3)
        u4 = unsampleline(l,p4)
        assert close(u0,-0.5)
        assert close(u1,0.0)
        assert close(u2,0.5)
        assert close(u3,1.0)
        assert close(u4,1.5)
        p = point(100,100)
        assert not unsampleline(l,p)

    def test_intersect(self):
        print("---> line-line intersection tests")
        a = point(5,0)
        b = point(0,5)
        c = point(-3,-3)
        d = point(1,1)
        e = point(10,10)
        l1 = line(a,b)
        print("l1:",vstr(l1))
        l2 = line(c,d)
        l3 = line(c,e)
        print("l1:" + vstr(l1) + ", l2:" + vstr(l2) +", l3:" + vstr(l3))

        int0 = lineLineIntersectXY(l1,l1)
        int0u = lineLineIntersectXY(l1,l1,params=True)
        int1 = lineLineIntersectXY(l1,l2,False)
        int1u = lineLineIntersectXY(l1,l2,params=True)
        int2 = lineLineIntersectXY(l1,l2,True)
        int2u = lineLineIntersectXY(l1,l2,params=True)
        int3 = lineLineIntersectXY(l1,l3,True)
        int3u = lineLineIntersectXY(l1,l3,params=True)
        
        assert not int0
        assert not int0u
        assert vclose(int1,point(2.5,2.5))
        assert len(int1u) == 2
        assert close(int1u[0],0.5)
        assert close(int1u[1],1.375)
        assert not int2
        assert close(int2u[0],0.5)
        assert close(int2u[1],1.375)
        assert vclose(int3,point(2.5,2.5))
        assert close(int3u[0], 0.5)
        assert close(int3u[1], 0.4230769230769231)
        
class TestArc:
    def test_create(self):
        arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0,-1)]
        arc11=arc(point(2.5,2.5),2.5,90.0,270.0)
        assert isarc(arc1)
        assert isarc(arc11)
        assert arc1 == arc11
        arc12=arc(arc11)
        assert isarc(arc12)
        assert arc12==arc11
        arc11[1][0]=5.0         # modify arc11
        assert arc12[1][0]==2.5 # assert that this didn't change arc12
        with pytest.raises(ValueError):
            assert not arc(vect(0,0),-2)

    def test_sample(self):
        arc1=arc(point(0,0),1.0,90.0,270.0)
        print("sample tests")
        print("arc1: ",vstr(arc1))
        for i in range(5):
            u = i/4
            theta = 90+i*45.0
            p = samplearc(arc1,u)
            x = point(cos(theta*pi2/360.0),sin(theta*pi2/360.0))
            print("sample: ",i," theta: ",theta)
            print("  sampled point p: ",vstr(p))
            print("  reference point x: ",vstr(x))
            print("  assert that these are the same")
            assert vclose(p,x)
        print("unsample tests")
        for i in range(5):
            u = i/4
            theta = 90+i*45.0
            x = point(cos(theta*pi2/360.0),sin(theta*pi2/360.0))
            uu = unsamplearc(arc1,x)
            print("unsample: ",i," theta: ",theta)
            print("  test point x: ",vstr(x))
            print("  unsampled parameter uu: ",uu)
            print("  calculated parameter u: ",u)
            print("  assert that these are the same")
            assert close(u,uu)

    def test_intersect_line(self):
        print("---> line-arc intersection testing")
        arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0)]
        print("arc1: {}".format(vstr(arc1)))
        a = point(5,0)
        b = point(0,5)
        c = point(-3,-3)
        d = point(0.0)
        l1 = line(a,b)
        print("l1: ", vstr(l1))
        l2 = line(c,d)
        print("l2: ",vstr(l2))
        int4 = lineArcIntersectXY(l1,arc1,False)
        int4u = lineArcIntersectXY(l1,arc1,params=True)
        print("cross-verifying intersection points and samples")
        for i in range(2):
            u1 = int4u[0][i]
            u2 = int4u[1][i]
            p=int4[i]
            p1 = sampleline(l1,u1)
            p2 = samplearc(arc1,u2)
            uu1 = unsampleline(l1,p)
            uu2 = unsamplearc(arc1,p)
            print("i: ",i,"u1: ",u1,"u2: ", u2,"p: ",vstr(p),
                  "p1: ",vstr(p1),"p2: ",vstr(p2))
            print("uu1: ",uu1," uu2: ",uu2)
            assert close(uu1,u1)
            assert close(uu2,u2)
            assert vclose(p,p1)
            assert vclose(p,p2)

        int5 = lineArcIntersectXY(l1,arc1,True)
        int6 = lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True)
        int7 = lineArcIntersectXY(l2,arc1,True)
        int8 = lineArcIntersectXY(l2,arc1,False)

        p1 = point(4.267766952966369, 0.7322330470336311)
        p2 = point(0.7322330470336311, 4.267766952966369)
        print("lineArcIntersectXY(l1,arc1,False): ", vstr(int4))
        assert len(int4) == 2
        assert vclose(int4[0],p1)
        assert vclose(int4[1],p2)
        print("lineArcIntersectXY(l1,arc1,params=True)",vstr(int4u))
        print("lineArcIntersectXY(l1,arc1,True): ", vstr(int5))
        assert len(int5) == 1
        assert vclose(int5[0],p2)
        print("tangent line, one intersection test")
        assert len(int6) == 1
        print("lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True): ",
              vstr(int6))
        assert vclose(int6[0],point(2.5,5.0))
        print("lineArcIntersectXY(l2,arc1,False): ",vstr(int7))
        assert not int7
        print("lineArcIntersectXY(l2,arc1,True): ",vstr(int8))
        assert len(int8) == 2
        assert vclose(int8[0], point(0.7322330470336311, 0.7322330470336311))
        assert vclose(int8[1], point(4.267766952966369, 4.267766952966369))


    def test_intersect_arc(self):
        print("---> arc-arc intersection tests")
        arc1=arc(vect(-2.5,0),2.5,90.0,270.0)
        arc2=arc(vect(2.5,0),2.5,90.0,270.0)

        print("arc1: {}".format(vstr(arc1)))
        print("arc2: {}".format(vstr(arc2)))
        print("single point of intersection test")
        int9 = arcArcIntersectXY(arc1,arc2,False)
        int9u =  arcArcIntersectXY(arc1,arc2,params=True)
        print("arcArcIntersectXY(arc1,arc2,False):",vstr(int9))
        print("arcArcIntersectXY(arc1,arc2,params=True):",int9u)
        
        assert len(int9) == 1
        assert vclose(int9[0],point(0,0))
        assert len(int9u) == 2
        assert close(int9u[0][0],-0.5)
        assert close(int9u[1][0],0.5)
        int10 = arcArcIntersectXY(arc1,arc2,True)
        int11 = arcArcIntersectXY(arc2,arc1,True)
        
        print("arcArcIntersectXY(arc1,arc2,True):",vstr(int10))
        print("arcArcIntersectXY(arc2,arc1,True):",vstr(int11))
        assert not int10
        assert not int11

        x = cos(pi/4)
        arc3=arc(vect(-x,-x),2.0,315.0,135.0)
        arc4=arc(vect(x,x),2.0,135.0,315.0)
        print("arc3: {}".format(vstr(arc3)))
        print("arc4: {}".format(vstr(arc4)))
        int12 = arcArcIntersectXY(arc3,arc4,True)
        int12u = arcArcIntersectXY(arc3,arc4,params=True)
        int13 = arcArcIntersectXY(arc3,arc4,False)
        print("arcArcIntersectXY(arc3,arc4,True):",vstr(int12))
        print("arcArcIntersectXY(arc3,arc4,False):",vstr(int13))
        print("arcArcIntersectXY(arc3,arc4,params=True):",int12u)
        assert len(int12u) == 2
        assert len(int12u[0])==2
        assert len(int12u[1])==2
        print("cross-verifying intersection points and samples")
        for i in range(2):
            u1 = int12u[0][i]
            u2 = int12u[1][i]
            p=int12[i]
            p1 = samplearc(arc3,u1)
            p2 = samplearc(arc4,u2)
            uu1 = unsamplearc(arc3,p)
            uu2 = unsamplearc(arc4,p)
            print("i: ",i,"u1: ",u1,"u2: ", u2,"p: ",vstr(p),
                  "p1: ",vstr(p1),"p2: ",vstr(p2))
            print("uu1: ",uu1," uu2: ",uu2)
            assert close(uu1,u1)
            assert close(uu2,u2)
            assert vclose(p,p1)
            assert vclose(p,p2)
            
        
class TestUtility:
    def test_copy(self):
        a = point(-5,-1)
        b = point(5,3)
        c = point(-3,-3)
        d = point(1,1)
        e = point(10,10)
        l = line(a,b)
        l2 = [c,d]
        l3 = line(c,e)
        foo = [a,b,l3,l2,d]
        foobar = deepcopy(foo)
        assert foobar == foo
        foo[0]=e
        assert not foobar==foo
        bar = [a,b,[1,2],l2,l3]
        
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

    def test_sample(self):
        a = point(0,-5)
        b = point(5,0)
        c = point(0,5)
        d = point(-5,0)

        pol1 = poly(a,b,c,d)
        print("---> tests for sampling and unsampling of polylines")
        print("pol1: ",vstr(pol1))
        pp = []
        print("sampled points: ")
        for i in range(9):
            u = i/8
            p = samplepoly(pol1,u)
            pp = pp+ [ p ]
            print("i: ",i,"u: ",u," p: ",vstr(p))
        print("unsampling polygon")
        for i in range(len(pp)):
            u = i/(len(pp)-1)
            uu = unsamplepoly(pol1,pp[i])
            print("i: ",i," p: ",vstr(pp[i])," u: ",u," uu: ",uu)
            assert close(u,uu)
            
    def test_intersect(self):
        a = point(0,-5)
        b = point(5,0)
        c = point(0,5)
        d = point(-5,0)
        
        pol1 = poly(a,b,c,d)
        pol2 = poly(a,b,c,d,a)

        line1 = line(point(0,0),point(5,5))
        line2 = line(point(0,-10),point(0,10))

        arc1 = arc(point(0,0),4.0,270,90)
        arc2 = arc(point(0,0),5.0,0,360)

        print("--> polyline and polygon intersection testing")
        print("pol1: ",vstr(pol1))
        print("pol2: ",vstr(pol2))
        print("line1: ",vstr(line1))
        print("line2: ",vstr(line2))
        print("arc1: ",vstr(arc1))
        print("arc2: ",vstr(arc2))
        
        int0 = intersectSimplePolyXY(line1,pol1,True)
        int1 = intersectSimplePolyXY(line1,pol1,False)
        int0u = intersectSimplePolyXY(line1,pol1,params=True)
        print("intersectSimplePolyXY(line1,pol1,True): ",vstr(int0))
        print("intersectSimplePolyXY(line1,pol1,False): ",vstr(int1))
        print("intersectSimplePolyXY(line1,pol1,params=True): ",vstr(int0u))
        assert len(int0) == 1
        assert vclose(int0[0],point(2.5,2.5))
        assert len(int1) == 1
        assert vclose(int1[0],point(2.5,2.5))
        assert len(int0u[0]) == 1
        assert close(int0u[0][0],0.5)
        assert close(int0u[1][0],0.5)

        int0 = intersectSimplePolyXY(line1,pol2,True)
        int1 = intersectSimplePolyXY(line1,pol2,False)
        int0u = intersectSimplePolyXY(line1,pol2,params=True)
        print("intersectSimplePolyXY(line1,pol2,True): ",vstr(int0))
        print("intersectSimplePolyXY(line1,pol2,False): ",vstr(int1))
        print("intersectSimplePolyXY(line1,pol2,params=True): ",vstr(int0u))
        assert len(int0) == 1
        assert vclose(int0[0],point(2.5,2.5))
        assert len(int1) == 2
        assert vclose(int1[0],point(2.5,2.5))
        assert vclose(int1[1],point(-2.5,-2.5))
        assert len(int0u[0]) == 2
        assert close(int0u[0][0],0.5)
        assert close(int0u[0][1],-0.5)
        assert close(int0u[1][0],0.375)
        assert close(int0u[1][1],0.875)

        int0 = intersectSimplePolyXY(arc1,pol2,True)
        int1 = intersectSimplePolyXY(arc1,pol2,False)
        int0u = intersectSimplePolyXY(arc1,pol2,params=True)
        print("intersectSimplePolyXY(arc1,pol2,True): ",vstr(int0))
        print("intersectSimplePolyXY(arc1,pol2,False): ",vstr(int1))
        print("intersectSimplePolyXY(arc1,pol2,params=True): ",vstr(int0u))


