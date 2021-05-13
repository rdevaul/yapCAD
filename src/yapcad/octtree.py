## quadtree and octtree representations for yapCAD geometry
## Born on 15 December 2020
## Copyright (c) 2020 Richard DeVaul

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""quadtree and octtree representations for yapCAD geometry"""

from yapcad.geom import *

# determine if two bounding boxes overlap.  This is a more challenging
# problem than it may appear at first glance.  It is not enough to
# test the corner points of one box to see if they fall inside the
# other, as illustrated by the following cases:

# case 1:  +--------------+
# inside   |     box1     |
#          |  +--------+  |
#          |  |  box2  |  |
#          |  +--------+  |
#          +--------------+

# No corner points of box1 lie inside box2, no lines intersect.

# case 2:     +------+
# cross       | box1 |
#     +-------+------+------+
#     | box2  |      |      |
#     +-------+------+------+
#             |      |
#             +------+

# No corner point of either box lie inside the other, projected lines
# intersect.

def boxoverlap2(bbx1,bbx2,dim3=True):
    """Determine if two bounding boxes overlap.  if dim3==True, treat the
    bounding boxes as 3D, otherwise treat them as co-planar 2D boxes.

    First, check to see if the maximum coordinates of one box are
    smaller than the minimum coordinates of the other, or vice versa.
    If so, no overlap is possible; return False

    if overlap is possible by test #1, check for the box-in-box
    special case for each box.  If so, return True

    Finally, for the 2D case: determine if horizontal lines of box1
    intersect with vertical lines of box2, and vice versa.  If any
    intersections found, return True, else return False

    For the 3D case, project the boxes into the XY, YZ, and XZ planes,
    and perform the 2D lines intersection check, as above.  Return
    True if and only if intersections are reported for each
    projection, otherwise return False
    """

    # check minmax
    if ((bbx1[1][0] < bbx2[0][0] and
         bbx1[1][1] < bbx2[0][1] and
         (not dim3 or bbx1[1][2] < bbx2[0][2])) or
        (bbx2[1][0] < bbx1[0][0] and
         bbx2[1][1] < bbx1[0][1] and
         (not dim3 or bbx2[1][2] < bbx1[0][2]))):
        return False # no overlap possible

    # check for box-in-box
    if ((bbx1[0][0] >= bbx2[0][0] and bbx1[1][0] <= bbx2[1][0] and
         bbx1[0][1] >= bbx2[0][1] and bbx1[1][1] <= bbx2[1][1] and
         (not dim3 or (bbx1[0][2] >= bbx2[0][2]
                       and bbx1[1][2] <= bbx2[1][2]))) or
        (bbx2[0][0] >= bbx1[0][0] and bbx2[1][0] <= bbx1[1][0] and
         bbx2[0][1] >= bbx1[0][1] and bbx2[1][1] <= bbx1[1][1] and
         (not dim3 or (bbx2[0][2] >= bbx1[0][2]
                       and bbx2[1][2] <= bbx1[1][2])))):
        return True

    def int2D(bb1,bb2,plane='XY'):
        """utility function for 2D box line intersection finding"""
        i = 0
        j = 0
        if plane == 'XY':
            i = 0
            j = 1
        elif plane == 'YZ':
            i = 1
            j = 2
        elif plane == 'XZ':
            i = 0
            j = 2
        else:
            raise ValueError('bad plane in int2D')

        def intHV(hln,vln):
            """ do a horizontal and vertial line intersect? """
            minx=hln[0][0]
            maxx=hln[1][0]
            if minx>maxx:
                swap = minx ; minx = maxx; maxx = swap

            if vln[0][0] < minx or vln[0][0] > maxx:
                return False
            miny=vln[0][1]
            maxy=vln[1][1]
            if miny > maxy:
                swap = miny ; miny = maxy ; maxy = swap
            if miny > hln[0][1] or maxy < hln[0][1]:
                return False
            return True
        
        # check for projected box-in-box
        if ((bb1[0][i] >= bb2[0][i] and bb1[1][i] <= bb2[1][i] and
             bb1[0][j] >= bb2[0][j] and bb1[1][j] <= bb2[1][j])
            or
            (bb2[0][i] >= bb1[0][i] and bb2[1][i] <= bb1[1][i] and
             bb2[0][j] >= bb1[0][j] and bb2[1][j] <= bb1[1][j])):
            return True

        # check for projected line intersections
        
        len1 = bb1[1][i] - bb1[0][i] # length
        wid1 = bb1[1][j] - bb1[0][j] # width
        len2 = bb2[1][i] - bb2[0][i] # length
        wid2 = bb2[1][j] - bb2[0][j] # width

        p0 = point(bb1[0][i],bb1[0][j])
        p1 = add(p0,point(len1,0))
        p2 = add(p1,point(0,wid1))
        p3 = add(p0,point(0,wid1))

        p4 = point(bb2[0][i],bb2[0][j])
        p5 = add(p4,point(len2,0))
        p6 = add(p5,point(0,wid2))
        p7 = add(p4,point(0,wid2))

        box1 = [[p0,p1],
                [p1,p2],
                [p2,p3],
                [p3,p0]]
        box2 = [[p4,p5],
                [p5,p6],
                [p6,p7],
                [p7,p4]]
        if intHV(box1[0],box2[1]):
            return True
        if intHV(box1[0],box2[3]):
            return True
        if intHV(box1[2],box2[1]):
            return True
        if intHV(box1[2],box2[3]):
            return True

        if intHV(box2[0],box1[1]):
            return True
        if intHV(box2[0],box1[3]):
            return True
        if intHV(box2[2],box1[1]):
            return True
        if intHV(box2[2],box1[3]):
            return True
        # if lineLineIntersectXY(box1[0],box2[1]):
        #     return True
        # if lineLineIntersectXY(box1[0],box2[3]):
        #     return True
        # if lineLineIntersectXY(box1[2],box2[1]):
        #     return True
        # if lineLineIntersectXY(box1[2],box2[3]):
        #     return True

        # if lineLineIntersectXY(box2[0],box1[1]):
        #     return True
        # if lineLineIntersectXY(box2[0],box1[3]):
        #     return True
        # if lineLineIntersectXY(box2[2],box1[1]):
        #     return True
        # if lineLineIntersectXY(box2[2],box1[3]):
        #     return True
        
        # for l1 in box1:
        #     for l2 in box2:
        #         if lineLineIntersectXY(l1,l2):
        #             return True
        return False
        
    # do projeted box line intersection tests
    return (int2D(bbx1,bbx2,'XY') and
            (not dim3 or int2D(bbx1,bbx2,'YZ')) and
            (not dim3 or int2D(bbx1,bbx2,'XZ')))
    
    
      
def boxoverlap(bbx1,bbx2,dim3=True):

    """determine if two bounding boxes overlap"""
    minx = bbx1[0][0]
    miny = bbx1[0][1]
    minz = bbx1[0][2]
    
    len1 = bbx1[1][0] - bbx1[0][0] # length
    wid1 = bbx1[1][1] - bbx1[0][1] # width
    hei1 = bbx1[1][2] - bbx1[0][2] # height

    p0=bbx1[0]
    p1=add(p0,point(len1,0,0))
    p2=add(p1,point(0,wid1,0))
    p3=add(p0,point(0,wid1,0))

    points = [p0,p1,p2,p3]
    
    if not dim3:
        for p in points:
            if isinsidebbox2D(bbx2,p):
                return True
        return False

    p4=add(p0,point(0,0,hei1))
    p5=add(p4,point(len1,0,0))
    p6=bbx1[1]
    p7=add(p4,point(0,wid1,0))

    points += [p4,p5,p6,p7]
    #print("boxoverlap points; ",vstr(points))
    #print("boxoverlap bbx2: ",vstr(bbx2))
    for p in points:
        if isinsidebbox(bbx2,p):
            return True
    return False
    
def bbox2oct(bbx,refbox,center):
    """
    Utility function to take a bounding box representation and assign
    it to zero or more octants.

    box2oct(bbx,refbox,center)

    bbx: 3D bounding box to assign
    refbox: reference 3D bounding box
    center: center point for purposes of assignment
      
    returns (potentially empty) list of octants, numbered 0 to 7
    """
    # print("bbox2oct :: bbx: ",vstr(bbx),"  refbox: ",vstr(refbox),"  center: ",vstr(center))
    rlist = []
    if not boxoverlap2(refbox,bbx):
        return rlist           # no overlap

    if bbx[0][2] < center[2]:
        rlist += bbox2quad(bbx,refbox,center)

    if bbx[1][2] >= center[2]:
        rlist += list(map(lambda x: x+4, bbox2quad(bbx,refbox,center)))

    return rlist
    
def bbox2quad(bbx,refbox,center):
    """Utility Function to take a bounding box representation and assign
       it to zero or more quads

    box2quad(bbx,refbox,center)

    bbx: 2D bounding box to assign
    refbox: reference 2D bounding box
    center: center point for purposes of assignment

    returns (potentially empty) list of quadrants, numbered 0 to 3
    """
    rlist = []
    if not boxoverlap2(refbox,bbx,dim3=False):
        return rlist           # no overlap

    if bbx[0][0] < center[0]:
        if bbx[0][1] < center[1]:
            rlist.append(2)
        if bbx[1][1] >= center[1]:
            rlist.append(1)
    if bbx[1][0] >= center[0]:
        if bbx[0][1] < center[1]:
            rlist.append(3)
        if bbx[1][1] >= center[1]:
            rlist.append(0)
    return rlist
            

def box2boxes(bbox,center,n,type='centersplit',elm=[]):
    """Function to take a bounding box (2d or 3D) and return a quad- or
    octtree decomposition of the box based on the value of center
    point, the type of split, and (potentially) the list of bounding
    boxes to be divvied up.
    """
    def boxmid(box):
        p0 = box[0]
        p1 = box[1]
        return scale3(add(p0,p1),0.5)

    # print("box2boxes :: bbox: ",vstr(bbox),"  center: ",vstr(center))
    if not isinsidebbox(bbox,center):
        raise ValueError('center point does not lie inside the bounding box')
    if type != 'centersplit':
        raise NotImplementedError('we are only doing center splits for now')
    if not n in (4,8):
        raise ValueError('bad tree dimension')

    cx = center[0]
    cy = center[1]
    cz = center[2]
        
    maxx = bbox[1][0]
    maxy = bbox[1][1]
    maxz = bbox[1][2]
    minx = bbox[0][0]
    miny = bbox[0][1]
    minz = bbox[0][2]

    if (maxx-minx < epsilon  or
        maxy-miny < epsilon or
        (n == 8 and maxz-minz < epsilon)) :
        raise ValueError('zero dimension box')
        
    if n == 4:
        z0 = 0
        z1 = 0
    else:
        z0 = minz
        z1 = cz
        
    box1 = [point(cx,cy,z0),
            point(maxx,maxy,z1)]
    box2 = [point(minx,cy,z0),
            point(cx,maxy,z1)]
    box3 = [point(minx,miny,z0),
            point(cx,cy,z1)]
    box4 = [point(cx,miny,z0),
            point(maxx,cy,z1)]

    r1 = [ [box1, boxmid(box1)],
           [box2, boxmid(box2)],
           [box3, boxmid(box3)],
           [box4, boxmid(box4)] ]

    if n == 4:
        return r1
    else:
        z0 = cz
        z1 = maxz
        
    box5 = [point(cx,cy,z0),
            point(maxx,maxy,z1)]
    box6 = [point(minx,cy,z0),
            point(cx,maxy,z1)]
    box7 = [point(minx,miny,z0),
            point(cx,cy,z1)]
    box8 = [point(cx,miny,z0),
            point(maxx,cy,z1)]

    r2 = [ [box5, boxmid(box5)],
           [box6, boxmid(box6)],
           [box7, boxmid(box7)],
           [box8, boxmid(box8)] ]

    return r1 + r2

def bboxdim(box):
    """ return length, width, and height of a bounding box"""

    length = box[1][0] - box[0][0]
    width = box[1][1] - box[0][1]
    height = box[1][2] - box[0][2]
    return [length, width, height]

def untag(tglist):
    """remove the (optional) tags that can be associated with elements in
    a geometry list."""
    glist = []
    for e in tglist:
        if isinstance(e,tuple):
            glist += [ e[0] ]
        else:
            glist += [ e ]
    return glist

class NTree():

    """Generalized n-tree representation for yapCAD geometry"""

    def __init__(self,n=8,geom=None,center=None,
                 mindim=None,maxdepth=None):

        if not n in [4,8]:
            raise ValueError('only quad- or octtrees supported')

        self.__n = n
        self.__depth = 0
        if mindim:
            if isinstance(mindim,(int,float)) and mindim > epsilon:
                self.__mindim = mindim
            else:
                raise ValueError('bad mindim value: '+str(mindim))
        else:
            self.__mindim = 1.0 # minimum dimension for tree element

        if maxdepth:
            if isinstance(maxdepth,int) and maxdepth > 0:
                self.__maxdepth = maxdepth
            else:
                raise ValueError('bad max depth value: '+str(maxdepth))
        elif self.__n == 4:
            self.__maxdepth = 8
        else: # n=8
            self.__maxdepth = 7

        self.__geom=[]
        if geom:
            ugeom = untag(geom)
            if isgeomlist(ugeom):
                self.__bbox= bbox(ugeom)
                self.__geom= geom
            else:
                raise ValueError('geom must be valid geometry list, or None')

        self.__center = None
        if center:
            if not ispoint(center):
                raise ValueError('center must be a valid point')
            self.__center = center
        else:
            if self.__geom != []:
                self.updateCenter()
                
        self.__tree = []
        self.__update=True

    def __repr__(self):
        return 'NTree(n={},depth={},mindim={},\n  geom={},\n  tree={})'.format(self.__n,self.__depth,self.__mindim,vstr(self.__geom),vstr(self.__tree))

    def addElement(self,element,tag=None):
        """ add a geometry element to the collection, don't update
        the tree -- yet """
        if not isgeomlist(element):
            raise ValueError('bad element passed to addElement')
        if not tag:
            gl = [ element ]
        else:
            gl = [ (element,tag) ]
        self.__geom += gl
        
        self.__update=True
        

    def updateCenter(self,center=None):
        """specify or compute new geometric center (or split poit) for tree,
        flag tree for rebuilding.

        """
        if self.__geom == []: #nothing to do
            return

        if not center:
            self.__center = scale3(add(self.__bbox[0],
                                       self.__bbox[1]),0.5)
        else:
            if ispoint(center):
                self.__center = center
            else:
                raise ValueError('bad center point passed to updateCenter')
        self.__update = True
        
        
    def updateTree(self):
        """
        build the tree from the current contents of the self.__geom list
        """
        if not self.__update or self.__geom == []:
            # nothing to do
            return
        self.__update = False

        # striptags for bbox
        ugeom = untag(self.__geom)
        self.__bbox = bbox(ugeom)
        if not self.__center:
            self.updateCenter()

        bxlist = list(map(lambda x: bbox(x), ugeom))
        bxidxlist = []
        for i in range(len(bxlist)):
            bxidxlist.append([i,bxlist[i]])
            
        self.__elem_idx_bbox = bxidxlist
        
        if not self.__center:
            self.updateCenter()

        # recursively build the tree
        def recurse(box,center,elements,depth=0):
            # print("recurse :: box: ",vstr(box),"  center: ",vstr(center))
            bbdim = bboxdim(box)
            if depth > self.__depth:
                self.__depth = depth
                
            if elements==[] or elements==None:
                return []
            elif (len(elements) <= self.__n or
                  depth > self.__maxdepth or
                  bbdim[0] < self.__mindim or
                  bbdim[1] < self.__mindim or
                  (self.__n == 8 and bbdim[2] < self.__mindim)):

                return [ 'e', box, center ] + list(map(lambda x:
                                                        x[0], elements))
            else:
                boxlist = ['b'] + box2boxes(box,center,self.__n)

                # print("boxlist: ",vstr(boxlist))
                # print("elements: ",vstr(elements))
                for e in elements:
                    func = None
                    box = e[1]
                    if self.__n == 8:
                        func = bbox2oct
                    else:
                        func = bbox2quad
                    ind = func(box,box,center)
                    # print("ind: ",ind,"  box: ",box,"  box: ",box,"  center: ",center)
                    #print ("e[0]: ",e[0], "  e[1]: ",vstr(e[1]),"  ind: ",ind)
                    # print ("box: ",vstr(box))
                    for j in ind:
                        boxlist[j+1].append(e)
                        
                for i in range(1,len(boxlist)):
                    boxlist[i] = recurse(boxlist[i][0],boxlist[i][1],
                                         boxlist[i][2:],depth+1)
                return boxlist
                
        self.__tree = recurse(self.__bbox, self.__center,
                              bxidxlist)
        # print("bxidxlist: ",vstr(bxidxlist))
        # print("self.__tree",self.__tree)

        return

    @property
    def depth(self):
        return self.__depth

    @depth.setter
    def depth(self,n):
        raise ValueError("can't set tree depth")

    @property
    def maxdepth(self):
        return self.__maxdepth

    @maxdepth.setter
    def maxdepth(self,d):
        if not isinstance(d,int) or d < 1:
            raise ValueError('bad maxdepth value: '+str(d))
        else:
            self.__maxdepth = d
    
    @property
    def mindim(self):
        return self.__mindim

    @mindim.setter
    def mindim(self,d):
        if not isinstance(d,(float,int)) or d < epsilon:
            raise ValueError('bad mindim value: '+str(d))
        else:
            self.__mindim = d
    
    def getElements(self,bbox):
        """return a list of geometry elements with bounding boxes that
        overalp the provided bounding box, or the empty list if none.

        """
        if self.__update:
            self.updateTree()

        bxidxlist = self.__elem_idx_bbox

        self.mxd = 0
        
        def recurse(subtree,depth=0):
            indices = []
            if depth > self.mxd:
                self.mxd = depth

            if subtree == None or subtree == []:
                return []
            elif subtree[0] == 'e':
                if boxoverlap2(subtree[1],bbox,dim3 = (self.__n==8)):
                    for ind in subtree[3:]:
                        box = bxidxlist[ind][1]
                        if boxoverlap2(bbox,box,dim3 = (self.__n==8)): 
                            indices.append(ind)
            else:
                for boxlist in subtree[1:]:
                    if len(subtree) == 1:
                        print('subtree ',subtree)
                    if len(boxlist) > 0:
                        indices += recurse(boxlist,depth+1)
            return indices

        idx = set(recurse(self.__tree))
                
        #print("unique indices: ",idx)

        elements = list(map(lambda x: self.__geom[x], list(idx)))

        return elements
