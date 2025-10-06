import pytest
import os
import random
from yapcad.geom import *
from yapcad.octtree import *

# Control flag for visual tests - set via environment variable or directly
VISUALTEST = os.environ.get('VISUALTEST', 'false').lower() in ('true', '1', 'yes')
if VISUALTEST:
    from yapcad.pyglet_drawable import pygletDraw

def three2two(x):
    """take a 3D bounding box with non-zero z coordinates, return
    a 2D bounding box with zero z coordinates"""
    a = vect(x[0])
    b = vect(x[1])
    a[2] = 0
    b[2] = 0
    return [a,b]

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


class TestOctree:
    """ Test utility functions """

    ## create non-overlappig bounding boxes
    b1 = [point(-10,-5,-5),point(-5,5,5)]
    b2 = [point(-2.5,-2.5,-2.5),point(2.5,2.5,2.5)]
    b3 = [point(5,-5,-5),point(10,5,5)]
    
    b4 = [point(-2.5,-10,-5),point(2.5,-5,5)]
    b5 = [point(-2.5,5,-5),point(2.5,10,5)]

    b6 = [point(-3.5,-3.5,-10),point(3.5,3.5,-5)]
    b7 = [point(-3.5,-3.5,5),point(3.5,3.5,10)]

    # a box that contains b2
    b10=[point(-3,-3,-3),point(3,3,3)]

    # crossing boxes that penetrate b2 and others

    b20= [point(-15,-1,-1),point(15,1,1)]
    b21= [point(-1,-15,-1),point(1,15,1)]
    b22= [point(-1,-1,-15),point(1,1,15)]

    # all of the boxes
    
    glist = [b1,b2,b3,b4,b5,b6,b7,
             b10,
             b20,b21,b22]

    # no 3D overlap list

    gl_no_over = [b1,b2,b3,b4,b5,b6,b7]

    # 3D crossing list

    gl_cross = [b10,
                b20,b21,b22]

    # no 2D overlap list
    gl2D_no_over = [b1,b2,b3,b4,b5]
    gl2D_no_over = list(map(three2two,gl2D_no_over))

    # 2D crossing list
    gl2D_cross = [b10, b20, b21]
    gl2D_cross = list(map(three2two,gl2D_cross))

    @pytest.mark.visual
    def test_boxoverlap2(self):
        """ do bounding box overlap check testing """
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()
        dd.linecolor = 'white'
        dd.cameradist = 35
    
        dd.draw_bbox(self.b1,dim3=True)
        dd.draw_bbox(self.b2,dim3=True)
        dd.draw_bbox(self.b3,dim3=True)
        dd.draw_bbox(self.b4,dim3=True)
        dd.draw_bbox(self.b5,dim3=True)
        dd.draw_bbox(self.b6,dim3=True)
        dd.draw_bbox(self.b7,dim3=True)

        dd.linecolor = 'red'
        dd.draw_bbox(self.b10,dim3=True)
        dd.linecolor = 'yellow'
        dd.draw_bbox(self.b20,dim3=True)
        dd.draw_bbox(self.b21,dim3=True)
        dd.draw_bbox(self.b22,dim3=True)
        
        dd.display()
    
        assert not boxoverlap2(self.b1,self.b2)
        assert not boxoverlap2(self.b1,self.b2,dim3=False)
        assert not boxoverlap2(self.b2,self.b3)
        assert not boxoverlap2(self.b2,self.b3,dim3=False)
        assert not boxoverlap2(self.b3,self.b4)
        assert not boxoverlap2(self.b3,self.b4,dim3=False)
        assert not boxoverlap2(self.b4,self.b2)
        assert not boxoverlap2(self.b4,self.b2,dim3=False)
        assert not boxoverlap2(self.b5,self.b2)
        assert not boxoverlap2(self.b5,self.b2,dim3=False)
    
        assert not boxoverlap2(self.b6,self.b2)
        # 2D box-in-box condition
        assert boxoverlap2(self.b6,self.b2,dim3=False)

        assert not boxoverlap2(self.b7,self.b2)
        # 2D box-in-box condition
        assert boxoverlap2(self.b7,self.b2,dim3=False)
        assert boxoverlap2(self.b2,self.b7,dim3=False)

        # 3D box-in-box condition
        assert boxoverlap2(self.b10,self.b2)
        assert boxoverlap2(self.b20,self.b2)
        assert boxoverlap2(self.b21,self.b2)
        assert boxoverlap2(self.b22,self.b2)

        bigbox = bbox(self.glist)
        # make sure computed geomlist bbox is correct, and assert that
        # all of the subboxes will be inside
        assert bigbox == [[-15,-15,-15,1],[15,15,15,1]]
        assert len(list(filter(lambda x: not boxoverlap2(bigbox,x),
                               self.glist))) == 0

    def test_bbox2oct(self):

        # find middle of box
        def boxmid(box):
            p0 = box[0]
            p1 = box[1]
            return scale3(add(p0,p1),0.5)

        bigbox = bbox(self.glist)
        center = boxmid(bigbox)

        inds = []
        for box in self.gl_no_over:
            inds.append(sorted(bbox2oct(box,bigbox,center)))

        # this is the correct list of octants for each test box
        assert inds[0] == [1,2,5,6]
        assert inds[1] == [0,1,2,3,4,5,6,7]
        assert inds[2] == [0,3,4,7]
        assert inds[3] == [2,3,6,7]
        assert inds[4] == [0,1,4,5]
        assert inds[5] == [0,1,2,3]
        assert inds[6] == [4,5,6,7]

        print("inds: ",inds)

    @pytest.mark.visual
    def test_quadtree(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()
        dd.cameradist = 35

        # make a quadtree
        bar = NTree(4)

        for l in self.gl2D_no_over:
            bar.addElement(l)

        bar.updateTree()
        dd.linecolor = 'red'
        for b in self.gl2D_no_over:
            dd.draw_bbox(b,dim3=False)
        dd.linecolor = 'white'
        dd.draw(self.gl2D_no_over)
        dd.display()

        # get full contents of tree
        bigbox = bbox(self.gl2D_no_over)
        contents = bar.getElements(bigbox)
        print("full contents: ",vstr(contents))
        assert len(contents) == len(self.gl2D_no_over)
        assert len(list(filter(lambda x: not x in self.gl2D_no_over,
                               contents))) == 0

        # retrieve the elements that overlap with self.b10

        contents = bar.getElements(self.b10)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 1
        cexpect = [three2two(self.b2)]
        for c in contents:
            assert c in cexpect

        # retrieve the elements that overlap with self.b20

        contents = bar.getElements(self.b20)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 3
        cexpect = list(map(three2two,[self.b1, self.b2, self.b3 ]))
        for c in contents:
            assert c in cexpect

               
        bar2 = NTree(4,mindim=4.0)
        # make 100 lines in space, keep track of which fall inside a
        # sub-region, then try to retrieve those with our tree,

        bigbox = [[-15,-15,0,1],[15,15,0,1]]
        deltabox = [[-10,-10,0,1],[10,10,0,1]]
        subbox = [[10,-10,0,1],[20,0,0,1]]
        glist=[]
        sublist=[]
        
        dim = 200
        plist1 = randomPoints(bigbox,dim)
        plist2 = randomPoints(deltabox,dim)
        
        for i in range(dim):
            l = line(plist1[i],add(plist1[i],plist2[i]))
            lbx = bbox(l)
            if boxoverlap2(lbx,subbox,dim3=False):
                sublist.append(l)
            else:
                glist.append(l)
            bar2.addElement(l)

        dd2 = pygletDraw()
        dd2.linecolor = [127,79,63]
        dd2.draw(glist)
        dd2.linecolor = 'white'
        dd2.draw_bbox(subbox,dim3=False)
        dd2.linecolor = 'aqua'
        dd2.draw(sublist)
        dd2.linecolor = 'red'
        for l in sublist:
            dd2.draw_bbox(l,dim3=False)
        dd2.display()

        bar2.updateTree()
        # print("bar2: ",bar2)
        # assert False

        print ("bar2 depth: ",bar2.depth)
        slist2 = bar2.getElements(subbox)
        assert len(slist2) == len(sublist)
        for l in slist2:
            assert l in sublist
 
    @pytest.mark.visual
    def test_octree(self):
        if not VISUALTEST:
            pytest.skip("Visual tests disabled (set VISUALTEST=true to enable)")

        dd = pygletDraw()
        dd.linecolor = 'white'
        dd.cameradist = 35

        # make an octree
        foo =NTree()
        # add the test boxes (lines) as elements

        for l in self.glist:
            foo.addElement(l)

#        for i in range(len(rp)):
#            foo.addElement(line(rp[i],rp2[i]))
    
        foo.updateTree()

        # draw the contents of the geometry list as lines
        dd.draw(self.glist)

        # draw the contents of the geometry list as bounding boxes

        dd.linecolor = 'red'
        for b in self.glist:
            dd.draw_bbox(b,dim3=True)
        dd.display()

        bigbox = bbox(self.glist)
        # get full contents of tree
        contents = foo.getElements(bigbox)
        print("full contents: ",vstr(contents))
        assert len(contents) == len(self.glist)
        assert len(list(filter(lambda x: not x in self.glist,
                               contents))) == 0

        # retrieve the elements that overlap with self.b10

        contents = foo.getElements(self.b10)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 5
        cexpect = [self.b2, self.b10, self.b20, self.b21, self.b22]
        for c in contents:
            assert c in cexpect

        # retrieve the elements that overlap with self.b1

        contents = foo.getElements(self.b1)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 2
        cexpect = [self.b1, self.b20 ]
        for c in contents:
            assert c in cexpect

        # retrieve the elements that overlap with self.b4

        contents = foo.getElements(self.b4)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 2
        cexpect = [self.b4, self.b21 ]
        for c in contents:
            assert c in cexpect

        # retrieve the elements that overlap with self.b7

        contents = foo.getElements(self.b7)
        print("selected contents: ",vstr(contents))

        assert len(contents) == 2
        cexpect = [self.b7, self.b22 ]
        for c in contents:
            assert c in cexpect
            
        
        foo2 = NTree(mindim=4.0)
        # make 1000 lines in space, keep track of which fall inside a
        # sub-region, then try to retrieve those with our tree,

        bigbox = [[-15,-15,-15,1],[15,15,15,1]]
        deltabox = [[-10,-10,-10,1],[10,10,10,1]]
        subbox = [[10,-10,5,1],[20,0,15,1]]
        glist=[]
        sublist=[]
        
        dim = 500
        plist1 = randomPoints(bigbox,dim)
        plist2 = randomPoints(deltabox,dim)
        
        for i in range(dim):
            l = line(plist1[i],add(plist1[i],plist2[i]))
            lbx = bbox(l)
            if boxoverlap2(lbx,subbox):
                sublist.append(l)
            else:
                glist.append(l)
            foo2.addElement(l)

        dd2 = pygletDraw()
        dd2.linecolor = [127,79,63]
        dd2.draw(glist)
        dd2.linecolor = 'white'
        dd2.draw_bbox(subbox,dim3=True)
        dd2.linecolor = 'aqua'
        dd2.draw(sublist)
        dd2.linecolor = 'red'
        for l in sublist:
            dd2.draw_bbox(l,dim3=True)
        dd2.display()

        foo2.updateTree()
            
        print ("foo2 depth: ",foo2.depth)
        slist2 = foo2.getElements(subbox)
        assert len(slist2) == len(sublist)
        for l in slist2:
            assert l in sublist
        
  
