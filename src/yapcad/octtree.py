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
    if not (boxoverlap(refbox,bbx) or boxoverlap(bbx,refbox)):
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
    if not (boxoverlap(refbox,bbx,dim3=False) or
            boxoverlap(bbx,refbox,dim3=False)):
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

class NTree():

    """Generalized n-tree representation for yapCAD geometry"""

    def __init__(self,n=8,geom=None,center=None):

        if not n in [4,8]:
            raise ValueError('only quad- or octtrees supported')

        self.__n = n

        self.__geom=[]
        if geom:
            if isgeomlist(geom):
                self.__bbox= bbox(geom)
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

    def addElement(self,element):
        """ add a geometry element to the collection, don't update
        the tree -- yet """
        gl = [ element ]
        if not isgeomlist(gl):
            raise ValueError('bad element passed to addElement')
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

        self.__bbox = bbox(self.__geom)
        if not self.__center:
            self.updateCenter()

        bxlist = list(map(lambda x: bbox(x), self.__geom))
        bxidxlist = []
        for i in range(len(bxlist)):
            bxidxlist.append([i,bxlist[i]])
            
        self.__elem_idx_bbox = bxidxlist
        
        if not self.__center:
            self.updateCenter()

        # recursively build the tree
        def recurse(bbox,center,elements):
            # print("recurse :: bbox: ",vstr(bbox),"  center: ",vstr(center))
            if elements==[]:
                return []
            elif len(elements) <= self.__n:
                return [ 'e', bbox, center ] + elements
            else:
                boxlist = ['b'] + box2boxes(bbox,center,self.__n)

                # print("boxlist: ",vstr(boxlist))
                # print("elements: ",vstr(elements))
                for e in elements:
                    func = None
                    if self.__n == 8:
                        func = bbox2oct
                    else:
                        func = bbox2quad
                    ind = func(e[1],bbox,center)
                    print ("e[0]: ",e[0], "  e[1]: ",vstr(e[1]),"  ind: ",ind)
                    # print ("bbox: ",vstr(bbox))
                    for i in ind:
                        boxlist[i+1].append(e[0])
                        
                for i in range(1,len(boxlist)):
                    boxlist[i] = recurse(boxlist[i][0],boxlist[i][1],
                                         boxlist[i][2:])
                return boxlist
                
        self.__tree = recurse(self.__bbox, self.__center,
                              bxidxlist)
        # print("bxidxlist: ",vstr(bxidxlist))
        print("self.__tree",self.__tree)
        self.__update = False

        return

    def getElements(self,bbox):
        """return a list of geometry elements with bounding boxes that
        overalp the provided bounding box, or the empty list if none.

        """
        if self.__update:
            self.updateTree()

        def recurse(subtree):
            flag = subtree[0]
            indices = []
            if subtree == []:
                pass
            elif flag == 'e':
                for elem in subtree[1:]:
                    if (boxoverlap(elem[1],bbox) or
                        boxoverlap(bbox,elem[1])):
                        indices.append(elem[0])
            else:
                for boxlist in subtree[1:]:
                    indices += recurse(boxlist)
            return indices

        idx = set(recurse(self.__tree))

        print("unique indices: ",idx)

        elements = list(map(lambda x: self.__geom[x], list(idx)))

        return elements
