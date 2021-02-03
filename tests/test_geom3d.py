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

    #vertices of tetrahedra
    tet1 = [ point(1, 1, 1),
             point(-1, -1, 1),
             point(-1, 1, -1),
             point(1, -1, -1) ]

    #alternate tetrahedron
    tet2 = [ point(-1, 1, 1),
             point(1, -1, 1),
             point(1, 1, -1),
             point(-1, -1, -1) ]

    #faces of tetrahedron 1
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

    #faces of tetrahedron 2
    face21 = [tet2[0],
              tet2[2],
              tet2[1]]
    face22 = [tet2[0],
              tet2[1],
              tet2[3]]
    face23 = [tet2[1],
              tet2[2],
              tet2[3]]
    face24 = [tet2[2],
              tet2[0],
              tet2[3]]
  
    def test_surface(self):
        """this test creates different types of closed two dimensional
        figures in the XY plane and turns them into surfaces by means
        of algorithmic tessellation.

        Tessellation is tricky, especially when dealing with
        non-convex shapes.  In this batch of tests we create
        non-convex figures using a variety of different yapCAD
        approaches

        The hope is that as tessellation bugs are found and fixed, new
        test cases will be added here to exercise those particular
        failure conditions, such that any regressions or new failure
        modes will be obvious.

        Holes in surfaces are currently not supported.
        """
        
        dd = pygletDraw()
        
        # The first test creates a sort of nautalus-shell-type shape
        # using a spiral.  We are starting with a representation using
        # arc segments, converting to a sampled polygon, and then
        # going from there to a surface.
        
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(-10,0),10.0,1.0)
        # "close" the arc with a line segment
        gl.append(line(sample(gl,1.0),
                       sample(gl,0.0)))

        dd.linecolor='silver'
        dd.draw(gl) #draw the original representation
        # turn the geometry list into a poly
        ply = geomlist2poly(gl)
        ply = translate(ply,point(0,0,0.5)) # pull it forward so that
                                            # when we draw it, it will
                                            # float slightly in front
                                            # of the rendered surface
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply) #draw the poly
        
        surf,bnd = poly2surfaceXY(ply) #make the surface
        dd.draw_surface(surf) # draw the surface

        # The second test is a big square non-convex "C" shape,
        # which should present a relatively simple test for the
        # surface algorithm
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
        dd.draw(ply2) # draw the poly
        surf2,bnd = poly2surfaceXY(ply2)
        dd.draw_surface(surf2) #draw the surface

        # The third test creates a similarly non-convex shape
        # as the previous test, but does so using rounded rect
        # Polygon instances and a boolean difference operation
        rr1 = makeRoundRect(5,5,0.5,center=point(-15,20))
        rr2 = makeRoundRect(6,3,0.5,center=point(-12,20))
        gl3 = Boolean('difference',[rr1,rr2])
        ply3 = geomlist2poly(gl3.geom())
        dd.linecolor='aqua'
        dd.draw(gl3)
        ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        surf3,bnd = poly2surfaceXY(ply3)
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply3)
        dd.draw_surface(surf3)

        # Finally, we make a big circular "C" that is nearly a
        # complete disk.  The large interior hole presents special
        # challenges for our surface tessellation algorithm

        # make the arcs as a geometry list -- note that the inside arc
        # is constructed in the samplereverse direction
        gl4 = [arc(point(0,7),10.0,0.0,350.0),
               arc(point(0,7),5.0,0.0,350.0,samplereverse=True)]
        dd.linecolor='aqua'
        dd.draw(gl4) #draw the arcs

        # convert our geometry list into a smapled polyline
        ply4 =geomlist2poly(gl4,minang = 10 )
        # close it into a polygon by appending a copy of the first
        # point to the end
        ply4.append(ply4[0])

        ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply4)

        #convert it into a surface
        surf4,bnd = poly2surfaceXY(ply4,minlen=1.0)
        dd.draw_surface(surf4)

        # turn our surface into a lines representation so we can
        # visualize the choices of the tessellation algorithm
        lines = surf2lines(surf4)
        lines = translate(lines,point(0,0,0.5))
        dd.linecolor='yellow'
        dd.draw(lines)

        # as a sanity check, the tesellation algorithm will return the
        # remaining boundry that it couldn't convert to a surface,
        # which might be zero-length, it might be degenerate (a couple
        # of points), or it could be a significant portion of geometry
        # that isn't converted due to some kind of bug.  Visualize
        # that.
        bnd.append(bnd[0])
        bnd = translate(bnd,point(0,0,2.0))
        dd.linecolor='white'
        dd.draw(bnd)
        
        print(f'len(surf4): {len(surf4)}')
        print(f'len vertices: {len(surf4[1])}')
        print(f'len normals: {len(surf4[2])}')
        print(f'len faces: {len(surf4[3])}')
        print(f'len boundaries: {len(surf4[4])}')
        print(f'len holes: {len(surf4[5])}')

        # check that our issurface function will correctly recognize
        # this surface
        assert issurface(surf4)
        assert issurface(surf4,fast=False)
        dd.display()

    
    def testFace(self):

        """method to test some computational geometry operations on faces.

        We create a tetrahedron, scale it, make faces, then separate
        those faces slightly from each other.

        We then generate a bunch of random test points, and compute
        the signed distance to the closest face.  For each test point,
        we also determine whether that point lies inside the
        tetrahedron, meaning that the signed distance to all faces is
        negative.

        We visualize the points as blue markers inside the
        tetrahedron, and red markers outside, with the color hue and
        luminance reflecting the closest distnace to a face.

        As an additonal test, we convert each face of the tetrahedron
        to a planar basis using tri2p0n, and then project each test
        point down into each plane.  If the signed distance is
        positive (the point lies outside the tetrahedron) and the
        projected point falls within the face, draw a line from that
        point to it's projected point in the face.  The result should
        be a bunch of face-orthogonal lines connecting to the test
        points.

        """
        dd = pygletDraw()

        # specify the verticies and colors of the tetrahedron faces.
        # Note, these faces are three coordinates, and thus not closed
        # polys.
        faces = [ self.face11,self.face12,
                  self.face13,self.face14]
        colors = [ 'red','yellow','green','blue' ]

        # lists for transformed face representations
        facetlist = [ [],[],[],[] ] # trandformation plane representation
        face2list= [ [],[],[],[] ] # in-face-plane coordinates

        # step through faces
        for i in range(len(faces)):
            face =faces[i]
            color = colors[i]
            face = scale(face,10.0) # scale by 10

            # compute the normal (and center, though we don't use it
            # here)
            p0,n = tri2p0n(face)
            face = translate(face,scale3(n,0.1)) # separate slightly
                                                 # from other faces
            # close poly to make bounding polygon
            face.append(face[0])
            faces[i] = face

            # convert face into planar basis
            facetlist[i] = tri2p0n(face,basis=True)
            p02,n2,tm,tmr = tuple(facetlist[i])
            print(f"p02: {p02}, n2 {n2}")
            print(f"tm: {tm}")

            # by definition, the origin should lie at the center point
            # of the face
            nz = tm.mul(p02)
            assert vclose(nz,point(0,0,0))

            # transform the 0,0,0 point back, and it should be the same
            # as the point we started wtih
            nzz = tmr.mul(nz)
            assert vclose(nzz,p02)

            # each vertex of the face, when transformed into the
            # planar basis, should lie in the z=0 plane.  Also, we will
            # memorize the transformed x-y coordinates of the vertices
            # for later use.
            face2 = []
            for p in face:
                p2 = tm.mul(p)
                face2.append(p2)
                # each point should fall in transformed z=0 plane
                print(f"p: {p}, p2: {p2}")
                assert close(p2[2],0.0)
            face2list[i] = face2

            # draw the face
            dd.linecolor = color
            dd.draw(face)

        #make surfaces -- this tests the creation of surfaces not
        #located in the XY plane
        for face in faces:
            print(f"face: {face}")
            surf,bnd = poly2surface(face)
            dd.draw_surface(surf)

        # specify surrounding box to contain random points, number of
        # points, and then generate those points
        dim = 400
        box = line(point(-10,-10,-10),
                   point(10,10,10))
        plist1 = randomPoints(box,dim)

        # for every point...
        for p in plist1:
            inside = True
            dmin = 100
            # and for every face
            for i in range(len(faces)):
                f= faces[i]
                f2 = face2list[i]
                
                # get the planar basis transform for the face
                ft = facetlist[i]
                # transform the point into this bsis
                p2 = ft[2].mul(p)
                # project the transformed point into the z=0 plane
                p2[2] = 0.0
                
                # compute the signed distance from original point to face
                d = signedFaceDistance(p,f)
                # (and unsigned)
                dm = abs(d)
                
                if d > 0: # if the signed distance is positive, see if
                          # the projected point lies inside the bounds
                          # of the face itself
                    if isInsideTriangleXY(p2,f2):
                        # if it does, then draw a line connecting the original
                        # point and the projected point
                        p3 = ft[3].mul(p2)
                        v = max(round(((7.0-d)*255/7.0)),0)
                        dd.linecolor=[0,v,255]
                        dd.draw(line(p3,p))
                    # record that the origial, untransformed point
                    # doesn't lie inside the tetrahedron
                    inside = False
                if dm < dmin: # if the signed distance is the smallest
                              # we've seen yet, then this is the new
                              # minimum distance
                    dmin = dm
            #compute a color factor to use to make red-ish and
            #blue-ish pints, depending on whether the point lies
            #inside or outside the tetrahedron.
            v = max(round(((7.0-dmin)*255/7.0)),0)
            # print(f"v: {v}")
            if inside:
                dd.linecolor=[0,v,255]
                dd.pointstyle='xo'
            else:
                dd.linecolor=[255,v,0]
                dd.pointstyle='o'

            #draw that point
            dd.draw(p)

        dd.display()
                             
