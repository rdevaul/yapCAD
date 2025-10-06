import pytest
import os
import random
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geom3d import *
from yapcad.geom3d_util import *
from yapcad.poly import *
from yapcad.combine import *

# Control flag for visual tests - set via environment variable or directly
VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw

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
              tet2[1],
              tet2[2]]
    face22 = [tet2[0],
              tet2[3],
              tet2[1]]
    face23 = [tet2[1],
              tet2[3],
              tet2[2]]
    face24 = [tet2[2],
              tet2[3],
              tet2[0]]
  
    @pytest.mark.visual
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

        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

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
        rr1 = RoundRect(5,5,0.5,center=point(-15,20))
        rr2 = RoundRect(6,3,0.5,center=point(-12,20))
        gl3 = Boolean('difference',[rr1,rr2])
        ply3 = geomlist2poly(gl3.geom)
        dd.linecolor='aqua'
        dd.draw(gl3)
        ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        surf3,bnd = poly2surfaceXY(ply3)
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply3)
        dd.draw_surface(surf3)

        # Finally, we make a disk. Surfaces with holes present a special
        # challenge to surface tessellation.

        # make the arcs as a geometry list -- note that the inside arc
        # is constructed in the samplereverse direction
        gl4 = [arc(point(0,7),8.0)]
        gl5 = [arc(point(0,7),5.0,samplereverse=True)]

        dd.linecolor='aqua'
        dd.draw(gl4) #draw the arcs
        dd.draw(gl5)

        # convert our geometry list into a smapled polyline
        ply4 =geomlist2poly(gl4,minang = 10 )
        ply5 =geomlist2poly(gl5,minang=10)
        # close it into a polygon by appending a copy of the first
        # point to the end
        #ply4.append(ply4[0])

        ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply4)

        #convert it into a surface
        surf4,bnd = poly2surfaceXY(ply4,holepolys=[ply5],minlen=1.0)
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
        if bnd:
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
        #import pdb ; pdb.set_trace()
        assert issurface(surf4,fast=False)
        dd.display()

    def test_extrude_surface_with_hole(self):
        outer = [
            point(0, 0),
            point(6, 0),
            point(6, 6),
            point(0, 6),
            point(0, 0),
        ]

        hole = [
            point(2, 2),
            point(4, 2),
            point(4, 4),
            point(2, 4),
            point(2, 2),
        ]

        surf, _ = poly2surfaceXY(outer, holepolys=[hole])
        solid = extrude(surf, 2.0)

        assert issolid(solid)

        bottom = solid[1][0]
        side = solid[1][1]

        total_edges = len(bottom[4]) + sum(len(loop) for loop in bottom[5])
        assert len(side[3]) == total_edges * 2

        offset = len(bottom[1])
        used_indices = {idx for face in side[3] for idx in face}
        for hole_loop in bottom[5]:
            for idx in hole_loop:
                assert idx in used_indices
                assert idx + offset in used_indices

    
    @pytest.mark.visual
    def test_face(self):
        """method to test some computational geometry operations on faces.
        """
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

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
                             
    @pytest.mark.visual
    def test_face_intersect(self):
        """Create two intersecting tetrahedra, draw them, then compute the
        lines representing the planar intersections"""
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()

        # specify the verticies and colors of the tetrahedron faces.
        # Note, these faces are three coordinates, and thus not closed
        # polys.
        faces1 = [ self.face11,self.face12,
                   self.face13,self.face14]
        colors1 = [ 'red','yellow','green','blue' ]

        faces2 = [ self.face21,self.face22,
                   self.face23,self.face24]
        colors2 = [ 'silver','gray','maroon','olive' ]

        # lists for transformed face representations
        facetlist1 = [] # trandformation plane representation
        face2list1= [] # in-face-plane coordinates
        facetlist2 = [] # trandformation plane representation
        face2list2= [] # in-face-plane coordinates

        def makeListsAndDraw(faces,colors):
            facetlist = []
            face2list = []
            # step through faces for 1st tetrahedron
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
                facetlist.append(tri2p0n(face,basis=True))
                p02,n2,tm,tmr = tuple(facetlist[i])

                # each vertex of the face, when transformed into the
                # planar basis, should lie in the z=0 plane.  Also, we will
                # memorize the transformed x-y coordinates of the vertices
                # for later use.
                face2 = list(map(lambda x: tmr.mul(x),face))
                face2list.append(face2)

                # draw the face
                dd.linecolor = color
                dd.draw(face)
            return facetlist,face2list

        
        facetlist1,face2list1 = makeListsAndDraw(faces1,colors1)
        facetlist2,face2list2 = makeListsAndDraw(faces2,colors2)

        #make surfaces -- this tests the creation of surfaces not
        #located in the XY plane
        for face in faces1 + faces2:
            surf,bnd = poly2surface(face)
            dd.draw_surface(surf)

        dd.linecolor='white'

        for f1 in faces1:
            for f2 in faces2:
                result = triTriIntersect(f1,f2,inside=True)
                print(f"result: {result}")
                if result != False:
                    dd.draw(result)

        triHard = [ point(-5,-5,15),
                    point(5,-5,15),
                    point(0,7,15),
                    point(-5,-5,15)]

        triHarder = [ point(-3,-3,15),
                    point(3,-3,15),
                    point(0,5,15),
                    point(-3,-3,15)]

        t0=triHard
        t1=rotate(triHard,180)
        t2=translate(triHard,point(-2,0,0))
        t3=translate(triHard,point(2,0,0))
        t4=triHarder
        
        dd.linecolor='blue'

        dd.draw(translate([t0,t1,t2,t3,t4],point(0,0,-.1)))

        dd.linecolor='white'
        i1 = triTriIntersect(t0,t4)
        assert i1
        dd.draw(i1)

        gg = intersectXY(t0,t1)
        i2 = triTriIntersect(t0,t1)
        print(f"t0: {t0}\nt1: {t1}")
        print(f"gg: {gg}")
        assert i2
        dd.linecolor='aqua'
        dd.draw(i2)

        i3 = triTriIntersect(t2,t3)
        assert i3
        dd.linecolor='white'
        dd.draw(i3)
        
        dd.display()


    @pytest.mark.visual
    def test_solid(self):
        """ Create solids procedurally and visualize them """
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()


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

        sld = extrude(surf,2.0)
        dd.make_object('obj1',lighting=True,linecolor='white',
                       material='jade',
                       position=point(0,0,0))
    
        dd.draw(sld,name='obj1')
        
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

        sld2 = extrude(surf2,5.0)
        dd.make_object('obj2',lighting=True,linecolor='white',
                       material='pearl',
                       position=point(0,0,0))
    
        dd.draw(sld2,name='obj2')

        # The third test creates a similarly non-convex shape
        # as the previous test, but does so using rounded rect
        # Polygon instances and a boolean difference operation
        rr1 = RoundRect(5,5,0.5,center=point(-15,20))
        rr2 = RoundRect(6,3,0.5,center=point(-12,20))
        gl3 = Boolean('difference',[rr1,rr2])
        ply3 = geomlist2poly(gl3.geom)
        dd.linecolor='aqua'
        dd.draw(gl3)
        ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        surf3,bnd = poly2surfaceXY(ply3)
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply3)

        sld3 = extrude(surf3,8.0)
        dd.make_object('obj3',lighting=True,linecolor='white',
                       material='obsidian',
                       position=point(0,0,0))
    
        dd.draw(sld3,name='obj3')


        # Finally, we make a disk. Surfaces with holes present a special
        # challenge to surface tessellation.

        # make the arcs as a geometry list -- note that the inside arc
        # is constructed in the samplereverse direction
        gl4 = [arc(point(0,7),8.0)]
        gl5 = [arc(point(0,7),5.0,samplereverse=True)]

        dd.linecolor='aqua'
        dd.draw(gl4) #draw the arcs
        dd.draw(gl5)

        # convert our geometry list into a smapled polyline
        ply4 =geomlist2poly(gl4,minang = 10 )
        ply5 =geomlist2poly(gl5,minang=10)
        # close it into a polygon by appending a copy of the first
        # point to the end
        #ply4.append(ply4[0])

        ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply4)

        #convert it into a surface
        surf4,bnd = poly2surfaceXY(ply4,holepolys=[ply5],minlen=1.0)

        #extrude it
        sld4 = extrude(surf4,5.0)
        dd.make_object('obj4',lighting=True,linecolor='white',
                       material='gold',
                       position=point(0,0,0))
    
        dd.draw(sld4,name='obj4')


        sld5 = sphere(10,center=point(10,0,20))
        dd.make_object('obj5',lighting=True,material='emerald',
                       position=point(0,0,0))
        dd.draw(sld5,name='obj5')

        sld6 = conic(5,5,7,center=point(-10,0,20))
        dd.make_object('obj6',lighting=True,material='ruby')
        dd.draw(sld6,name='obj6')

        ## do a quick revolution surface test

        def fcontour(x,ir=10,len=20,fr=5):
            if x < 0:
                return 0
            elif x <= ir:
                xx = ir-x
                return sqrt(ir*ir-xx*xx)
            elif x <= len+ir:
                d = len+ir-x
                u = d/len
                return u*ir+(1-u)*fr
            elif x <= len+ir+fr:
                xx = x-(ir+len)
                return sqrt(fr*fr-xx*xx)
            else:
                return 0
            
        revs = makeRevolutionSurface(fcontour,0,35,50)
        revs2 = deepcopy(revs)
        revs = translatesurface(revs,point(0,0,30))
        revs = rotatesurface(revs,90,cent=point(0,0,40),axis=point(0,1,0))
        airship = solid([revs],[],['procedure',
                                   f"makeRevolutionSurface(fcontour,0,35,50)"])
        dd.make_object('airship',lighting=True,material='black rubber')
        dd.draw(airship,name='airship')

        airship2 = solid([revs2],
                         [],['procedure',
                             f"makeRevolutionSurface(fcontour,0,35,50)"])
        airship2 = rotate(airship2,90,cent=point(0,0,17.5),axis=vect(0,1,0,0))
        airship2 = mirror(airship2,"yz")
        airship2 = translate(airship2,point(0,0,50))
        dd.draw(airship2,name='airship')

        dd.display()


class TestSolidTopology:
    """Test solid topology analysis functions: issolidclosed and volumeof"""

    def test_issolidclosed_empty_solid(self):
        """Empty solid should be considered closed"""
        empty = solid([])
        assert issolidclosed(empty)

    def test_issolidclosed_single_surface(self):
        """Single open surface is not a closed solid"""
        # Create a simple triangular surface (not closed)
        p1 = point(0, 0, 0)
        p2 = point(1, 0, 0)
        p3 = point(0.5, 1, 0)
        surf = surface([p1, p2, p3],
                      [[0, 0, 1, 0]] * 3,
                      [[0, 1, 2]])

        s = solid([surf])
        # Single surface can't be closed - edges not shared by 2 faces
        assert not issolidclosed(s)

    def test_issolidclosed_tetrahedron(self):
        """Tetrahedron is a closed solid"""
        # Create tetrahedron vertices
        tet_verts = [
            point(1, 1, 1),
            point(-1, -1, 1),
            point(-1, 1, -1),
            point(1, -1, -1)
        ]

        # Create four triangular faces with correct winding
        faces = [
            [0, 2, 1],  # face pointing outward
            [0, 1, 3],
            [1, 2, 3],
            [2, 0, 3]
        ]

        # Compute normals for each vertex of each face
        normals = []
        for f in faces:
            tri = [tet_verts[f[0]], tet_verts[f[1]], tet_verts[f[2]]]
            p0, n = tri2p0n(tri)
            normals.extend([n, n, n])

        # Create vertices list - each vertex appears in multiple faces
        all_verts = []
        all_norms = []
        for f_idx, f in enumerate(faces):
            for v_idx in f:
                all_verts.append(tet_verts[v_idx])
                all_norms.append(normals[f_idx * 3])

        # Remap face indices to the expanded vertex list
        new_faces = []
        for i in range(len(faces)):
            new_faces.append([i * 3, i * 3 + 1, i * 3 + 2])

        surf = surface(all_verts, all_norms, new_faces)
        tet = solid([surf])

        assert issolidclosed(tet)

    def test_issolidclosed_cube(self):
        """Cube created by prism should be closed"""
        from yapcad.geom3d_util import prism
        cube = prism(2, 2, 2)
        assert issolidclosed(cube)

    def test_issolidclosed_sphere(self):
        """Sphere should be closed"""
        from yapcad.geom3d_util import sphere
        sph = sphere(2.0, depth=2)
        assert issolidclosed(sph)

    def test_issolidclosed_extruded_square(self):
        """Extruded square should be closed"""
        from yapcad.geom3d_util import extrude
        # Create square polygon
        square = [point(0, 0), point(2, 0), point(2, 2), point(0, 2), point(0, 0)]
        surf, _ = poly2surfaceXY(square)
        extruded = extrude(surf, 3.0)
        assert issolidclosed(extruded)

    def test_issolidclosed_extruded_with_hole(self):
        """Extruded surface with hole should be closed"""
        from yapcad.geom3d_util import extrude
        # Outer square
        outer = [point(0, 0), point(4, 0), point(4, 4), point(0, 4), point(0, 0)]
        # Inner square (hole)
        hole = [point(1, 1), point(3, 1), point(3, 3), point(1, 3), point(1, 1)]

        surf, _ = poly2surfaceXY(outer, holepolys=[hole])
        extruded = extrude(surf, 2.0)
        assert issolidclosed(extruded)

    def test_issolidclosed_open_box(self):
        """Box missing a face should not be closed"""
        from yapcad.geom3d_util import rectangularPlane, extrude

        # Create bottom surface
        bottom = rectangularPlane(2, 2, point(0, 0, 0))

        # Create four side surfaces manually (missing top)
        side1 = rectangularPlane(2, 1, point(0, 0, 0.5))
        side1 = rotatesurface(side1, 90, axis=point(1, 0, 0))
        side1 = translatesurface(side1, point(0, -1, 0))

        side2 = rectangularPlane(2, 1, point(0, 0, 0.5))
        side2 = rotatesurface(side2, 90, axis=point(1, 0, 0))
        side2 = translatesurface(side2, point(0, 1, 0))

        side3 = rectangularPlane(2, 1, point(0, 0, 0.5))
        side3 = rotatesurface(side3, 90, axis=point(0, 1, 0))
        side3 = translatesurface(side3, point(-1, 0, 0))

        side4 = rectangularPlane(2, 1, point(0, 0, 0.5))
        side4 = rotatesurface(side4, 90, axis=point(0, 1, 0))
        side4 = translatesurface(side4, point(1, 0, 0))

        # Create open box (no top)
        open_box = solid([bottom, side1, side2, side3, side4])

        # Should not be closed
        assert not issolidclosed(open_box)

    def test_volumeof_unit_cube(self):
        """Volume of unit cube should be 1"""
        from yapcad.geom3d_util import prism
        cube = prism(1, 1, 1)
        vol = volumeof(cube)
        assert abs(vol - 1.0) < 0.001

    def test_volumeof_cube_2x2x2(self):
        """Volume of 2x2x2 cube should be 8"""
        from yapcad.geom3d_util import prism
        cube = prism(2, 2, 2)
        vol = volumeof(cube)
        assert abs(vol - 8.0) < 0.001

    def test_volumeof_rectangular_prism(self):
        """Volume of 3x4x5 prism should be 60"""
        from yapcad.geom3d_util import prism
        box = prism(3, 4, 5)
        vol = volumeof(box)
        assert abs(vol - 60.0) < 0.01

    def test_volumeof_sphere(self):
        """Volume of sphere should match 4πr³/3"""
        from yapcad.geom3d_util import sphere
        import math

        radius = 2.0
        expected = (4.0 / 3.0) * math.pi * (radius ** 3)

        sph = sphere(2 * radius, depth=3)  # depth=3 for better approximation
        vol = volumeof(sph)

        # Allow 1% error due to triangulation approximation
        assert abs(vol - expected) / expected < 0.01

    def test_volumeof_extruded_square(self):
        """Volume of extruded 2x2 square by height 3 should be 12"""
        from yapcad.geom3d_util import extrude

        square = [point(0, 0), point(2, 0), point(2, 2), point(0, 2), point(0, 0)]
        surf, _ = poly2surfaceXY(square)
        extruded = extrude(surf, 3.0)

        vol = volumeof(extruded)
        assert abs(vol - 12.0) < 0.01

    def test_volumeof_extruded_with_hole(self):
        """Volume of extruded surface with hole"""
        from yapcad.geom3d_util import extrude

        # 4x4 outer square with 2x2 inner hole
        # Expected volume: (4*4 - 2*2) * height = 12 * 2 = 24
        outer = [point(0, 0), point(4, 0), point(4, 4), point(0, 4), point(0, 0)]
        hole = [point(1, 1), point(3, 1), point(3, 3), point(1, 3), point(1, 1)]

        surf, _ = poly2surfaceXY(outer, holepolys=[hole])
        extruded = extrude(surf, 2.0)

        vol = volumeof(extruded)
        expected = (4 * 4 - 2 * 2) * 2
        assert abs(vol - expected) < 0.1

    def test_volumeof_truncated_cone(self):
        """Volume of truncated cone (frustum) should match formula"""
        from yapcad.geom3d_util import conic
        import math

        base_radius = 3.0
        top_radius = 1.0
        height = 4.0
        # Volume of frustum: (π*h/3)(r1² + r1*r2 + r2²)
        expected = (math.pi * height / 3.0) * (
            base_radius**2 + base_radius*top_radius + top_radius**2
        )

        frustum = conic(base_radius, top_radius, height, angr=10)
        vol = volumeof(frustum)

        # Allow 2% error due to triangulation
        assert abs(vol - expected) / expected < 0.02

    def test_volumeof_cylinder(self):
        """Volume of cylinder should match πr²h"""
        from yapcad.geom3d_util import conic
        import math

        radius = 2.0
        height = 5.0
        expected = math.pi * (radius ** 2) * height

        cylinder = conic(radius, radius, height, angr=5)
        vol = volumeof(cylinder)

        # Allow 1% error
        assert abs(vol - expected) / expected < 0.01

    def test_volumeof_requires_closed_solid(self):
        """volumeof should raise error for non-closed solid"""
        # Create a single triangular surface (not closed)
        p1 = point(0, 0, 0)
        p2 = point(1, 0, 0)
        p3 = point(0.5, 1, 0)
        surf = surface([p1, p2, p3],
                      [[0, 0, 1, 0]] * 3,
                      [[0, 1, 2]])

        s = solid([surf])

        with pytest.raises(ValueError, match="must be topologically closed"):
            volumeof(s)

    def test_volumeof_tetrahedron(self):
        """Volume of tetrahedron should match a³/(6√2)"""
        import math

        # Create a tetrahedron with vertices at specific positions
        # This forms a tetrahedron with all edges of equal length
        tet_verts = [
            point(1, 1, 1),
            point(-1, -1, 1),
            point(-1, 1, -1),
            point(1, -1, -1)
        ]

        # Create faces
        faces = [
            [0, 2, 1],
            [0, 1, 3],
            [1, 2, 3],
            [2, 0, 3]
        ]

        # Create surface
        all_verts = []
        all_norms = []
        for f_idx, f in enumerate(faces):
            tri = [tet_verts[f[0]], tet_verts[f[1]], tet_verts[f[2]]]
            p0, n = tri2p0n(tri)
            for v_idx in f:
                all_verts.append(tet_verts[v_idx])
                all_norms.append(n)

        new_faces = []
        for i in range(len(faces)):
            new_faces.append([i * 3, i * 3 + 1, i * 3 + 2])

        surf = surface(all_verts, all_norms, new_faces)
        tet = solid([surf])

        # Verify it's closed
        assert issolidclosed(tet)

        # Calculate actual edge length using the length() function
        edge = line(tet_verts[0], tet_verts[1])
        a = length(edge)

        # Volume of regular tetrahedron: V = a³/(6√2)
        expected = (a ** 3) / (6 * math.sqrt(2))

        vol = volumeof(tet)

        # Check volume matches expected within small tolerance
        assert abs(vol - expected) < 0.001
