import pytest
import random
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geom3d import *
from yapcad.pyglet_drawable import *
from yapcad.poly import *
from yapcad.combine import *

def randomPoints(bbox,numpoints):
    """Given a 3D bounding box and a number of points to generate, 
    return a list of uniformly generated random points within the 
    bounding box"""
    
    points = []
    minx = bbox[0][0]
    maxx = bbox[1][0]
    miny = bbox[0][1]
    maxy = bbox[1][1]
    minz = bbox[0][2]
    maxz = bbox[1][2]
    rangex = maxx-minx
    rangey = maxy-miny
    rangez = maxz-minz
    for i in range(numpoints):
        points.append(point(random.random()*rangex+minx,
                            random.random()*rangey+miny,
                            random.random()*rangez+minz))
    return points


class TestSurface:
    """test surface representation"""
    def test_surface(self):
        
        dd = pygletDraw()
        dd.linecolor='silver'
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(-10,0),10.0,1.0)
        # "close" the arc with a line segment
        gl.append(line(sample(gl,1.0),
                       sample(gl,0.0)))
        dd.draw(gl)
        # turn the geometry list into a poly

        ply = geomlist2poly(gl)
        ply = translate(ply,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply)
        surf,bnd = poly2surface(ply)
        dd.draw_surface(surf)

        ply2 = [ point(20,20),
                 point(10,20),
                 point(10,-20),
                 point(20,-20),
                 point(20,-10),
                 point(15,-10),
                 point(15,10),
                 point(20,10),
                 point(20,20) ]
                 
        ply2 = translate(ply2,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply2)
        surf2,bnd = poly2surface(ply2)
        dd.draw_surface(surf2)

        rr1 = makeRoundRect(5,5,0.5,center=point(-15,20))
        rr2 = makeRoundRect(6,3,0.5,center=point(-12,20))
        gl3 = Boolean('difference',[rr1,rr2])
        ply3 = geomlist2poly(gl3.geom())
        dd.linecolor='aqua'
        dd.draw(gl3)
        ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        surf3,bnd = poly2surface(ply3)
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply3)
        dd.draw_surface(surf3)

        gl4 = [arc(point(0,7),10.0,0.0,270.0),
               arc(point(0,7),5.0,0.0,270.0,samplereverse=True)]
        gl4 = [arc(point(0,7),10.0,0.0,350.0),
               arc(point(0,7),5.0,0.0,350.0,samplereverse=True)]
        dd.linecolor='aqua'
        dd.draw(gl4)
        ply4 =geomlist2poly(gl4,minang = 10 )
        ply4.append(ply4[0])

        ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply4)

        surf4,bnd = poly2surface(ply4,minlen=1.0)
        # assert False
        lines = surf2lines(surf4)
        lines = translate(lines,point(0,0,0.5))
        dd.draw_surface(surf4)
        dd.linecolor='yellow'
        dd.draw(lines)
        bnd.append(bnd[0])
        bnd = translate(bnd,point(0,0,1.0))
        dd.linecolor='white'
        dd.draw(bnd)
        print(f'len(surf4): {len(surf4)}')
        print(f'len vertices: {len(surf4[1])}')
        print(f'len normals: {len(surf4[2])}')
        print(f'len faces: {len(surf4[3])}')
        print(f'len boundaries: {len(surf4[4])}')
        print(f'len holes: {len(surf4[5])}')
        assert issurface(surf4)
        assert issurface(surf4,fast=False)
        dd.display()

    
    def testFace(self):
        """
        method to test computational geometry operations on faces
        """
        tet1 = [ point(1, 1, 1),
                 point(-1, -1, 1),
                 point(-1, 1, -1),
                 point(1, -1, -1) ]

        tet2 = [ point(-1, 1, 1),
                 point(1, -1, 1),
                 point(1, 1, -1),
                 point(-1, -1, -1) ]
        
        dd = pygletDraw()
        dd.linecolor='silver'

        face11 = [tet1[0],
                  tet1[2],
                  tet1[1]]
        face12 = [tet1[0],
                  tet1[1],
                  tet1[3]]
        face13 = [tet1[1],
                  tet1[2],
                  tet1[3]]
        face14 = [tet1[2],
                  tet1[0],
                  tet1[3]]

        faces = [ face11,face12,face13,face14]
        colors = [ 'red','yellow','green','blue' ]
        facetlist = [ [],[],[],[] ]
        face2list= [ [],[],[],[] ]
        for i in range(len(faces)):
            face =faces[i]
            color = colors[i]
            p0,n = tuple(tri2p0n(face))
            face = scale(face,10.0)
            face = translate(face,scale3(n,0.1))
            face.append(face[0])
            dd.linecolor = color
            faces[i] = face

            facetlist[i] = tri2p0n(face,basis=True)
            p02,n2,tm,tmr = tuple(facetlist[i])
            print(f"p02: {p02}, n2 {n2}")
            print(f"tm: {tm}")
            nz = tm.mul(p02)
            assert vclose(nz,point(0,0,0))
            nzz = tmr.mul(nz)
            assert vclose(nzz,p02)
            face2 = []
            for p in face:
                p2 = tm.mul(p)
                face2.append(p2)
                # each point should fall in transformed z=0 plane
                print(f"p: {p}, p2: {p2}")
                assert close(p2[2],0.0)
            face2list[i] = face2
            dd.draw(face)

        dim = 400
        box = line(point(-10,-10,-10),
                   point(10,10,10))
        plist1 = randomPoints(box,dim)

        dd.linecolor = 'silver'
        for p in plist1:
            inside = True
            dmin = 100
            for i in range(len(faces)):
                f= faces[i]
                f2 = face2list[i]
                ft = facetlist[i]
                p2 = ft[2].mul(p)
                p2[2] = 0.0
                d = signedFaceDistance(p,f)
                dm = abs(d)
                if d > 0:
                    if isInsideTriangleXY(p2,f2):
                        p3 = ft[3].mul(p2)
                        dd.draw(line(p3,p))
                    inside = False
                if dm < dmin:
                    dmin = dm
            v = max(round(((7.0-dmin)*255/7.0)),0)
            # print(f"v: {v}")
            if inside:
                dd.linecolor=[0,v,255]
                dd.pointstyle='xo'
            else:
                dd.linecolor=[255,v,0]
                dd.pointstyle='o'
            dd.draw(p)

        dd.display()
                             
