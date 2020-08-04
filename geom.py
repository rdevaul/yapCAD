## simple computational geometry library for yapCAD
## Born on 29 July, 2020
## Richard DeVaul

from math import *

## constants
epsilon=0.0000001
pi2 = 2.0*pi

## operations on vectors
## ------------------------

## vectors are defined as a list of four numbers, i.e. [x,y,z,w]

## The use of four coordinates supports the generalized homogemeous
## coordinates or projective coordinates popular with computer
## graphicists in the 1990s. In general, the w coordinate is a
## normalization coordinate, and is either 1.0 or all coordinates are
## assumed to be scaled by w.  Simple vector operations are assumed to
## operate in the w=1 hyperplane, and no checking of w is performed.
## See https://en.wikipedia.org/wiki/Homogeneous_coordinates

## utility function to determine if argument is a "real" python
## number, since booleans are considered ints (True=1 and False=0 for
## integer arithmatic) but 1 and 0 are not considered boolean

def isgoodnum(n):
    return (not isinstance(n,bool)) and isinstance(n,(int,float))

## Convenience function for making a homogeneous coordinates 4 vector
## from practically anything

def vect(a=False,b=False,c=False):
    r = [0,0,0,1]
    if isgoodnum(a):
        r[0]=a
        if isgoodnum(b):
            r[1]=b
            if isgoodnum(c):
                r[2]=c
    elif isinstance(a,(tuple,list)):
        for i in range(min(3,len(a))):
            x=a[i]
            if isgoodnum(x):
                r[i]=x
    return r

## check to see if argument is a proper vector for our purposes
def isvect(x):
    return isinstance(x,list) and len(x) == 4 and isgoodnum(x[0]) and isgoodnum(x[1]) and isgoodnum(x[2]) and isgoodnum(x[3])
    
## R^3 -> R^3 functions: ignore w component
## ------------------------------------------------
def add(a,b):
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2],1.0]

def sub(a,b):
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2],1.0]

def scale(a,c):
    return [a[0]*c,a[1]*c,a[2]*c,1.0]

## NOTE: this function assumes that a lies in the x,y plane.  If this
## is not the case, the results are bogus.
def ortho(a):
    # compute an orthogonal vector to a, by
    # crossing [a1, a2, a3] with [0, 0, 1]

    return [ a[1], -a[0], 0, 1.0 ]

## Compute the cross generalized product of a x b, assuming that both
## fall into the w=1 hyperplane
def cross(a,b):

    return [ a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0],
             1.0 ]

## R^4 -> R^4 functions: operate on w component
def add4(a,b):
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2],a[3]+b[3]]

def sub4(a,b):
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2],a[3]-b[3]]

def scale4(a,c):
    return [a[0]*c,a[1]*c,a[2]*c,a[3]*c]

## Homogenize, or project back to the w=1 plane by scaling all values
## by w
def homo(a):
    return [ a[0]/a[3],
             a[1]/a[3],
             a[2]/a[3],
             1 ]

## R^3 -> R functions -- ignore w component
## ----------------------------------------
def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def mag(a):
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

def dist(a,b):  # compute distance between two points a & b
    return mag(sub(a,b))

## R^3 -> bool functions
## ---------------------

# utility function to determine if a list of points lies in one of the
# cardinal planes: XY, YZ, XZ

def isCardinalPlanar(plane="xy",points=[]):
    if plane=="xy":
        idx = 2
    elif plane == "yz":
        idx = 0
    elif plane == "xz":
        idx = 1
    else:
        raise ValueError('bad plane specification: {}'.format(plane))
        return False

    if len(points) < 2:
        return True

    ref = points[0][idx]
    for i in range(1,len(points)):
        if abs(points[i][idx]-ref) > epsilon:
            return False
    return True

## conveniience functions
def isXYPlanar(points=[]):
    return isCardinalPlanar("xy",points)

def isYZPlanar(points=[]):
    return isCardinalPlanar("yz",points)

def isXZPlanar(points=[]):
    return isCardinalPlanar("xz",points)


## R^4 -> R functions 
## ----------------------------------------
def dot4(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

def mag4(a):
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3])

def dist4(a,b):  # compute distance between two points a & b
    return mag(sub4(a,b))


## misc operations
## -----------------------------------------

# pretty printing string formatter for vectors, lines, and polygons.
# You can use this anywhere you use str(), since it will fall back to
# str() if the argument isn't a yapCAD vector, line, or polygon.
def vstr(a):
    # utility functions for recursively checking and formatting lists
    def _isallvect(foo):
        if len(foo)==1:
            return isvect(foo[0])
        else:
            return isvect(foo[0]) and _isallvect(foo[1:])
    def _makestr(foo):
        if len(foo) ==1:
            return vstr(foo[0])
        else:
            return vstr(foo[0]) + ", " + _makestr(foo[1:])
    # if it's not a list, it's not a vector, line, or poly    
    if not isinstance(a,list):
        return str(a) # fall back to default string formatting
    # if it is a vector, format it to leave out extraneous
    # coordinates.  NOTE: this means that 3 vectors that happen to
    # fall into the z=0 plane will be formatted as though they were 2
    # vectors.
    if isvect(a):
        if abs(a[3]-1.0) > epsilon: # not in w=1
            return "[{}, {}, {}, {}]".format(a[0],a[1],a[2],a[3])
        elif abs(a[2]) > epsilon: # not in z=0
            return "[{}, {}, {}]".format(a[0],a[1],a[2])
        else: # in x-y plane
            return "[{}, {}]".format(a[0],a[1])
    # if it is a list of vectors (line or poly) format it appropriately
    elif len(a)>0 and _isallvect(a):
        return "["+_makestr(a)+"]"
    else: #it's none of those things, fall back to default string
          #formatting
        return str(a)
    

## operations on lines
## --------------------
## lines are defined as lists of two points, i.e.
## [vect(x1, y1),vect(x2, y2)]

## Compute the intersection of two lines that lie in the same x,y plane
def intersectXY(l1,l2,inside=True):
    x1=l1[0][0]
    y1=l1[0][1]
    z1=l1[0][2]
    
    x2=l1[1][0]
    y2=l1[1][1]
    z2=l1[1][2]

    x3=l2[0][0]
    y3=l2[0][1]
    z3=l2[0][2]
    
    x4=l2[1][0]
    y4=l2[1][1]
    z4=l2[1][2]

    ## check for x,y planar consistency
    if abs(z2-z1) > epsilon or abs(z3-z1) > epsilon or abs(z4-z1) > epsilon:
        return False

    ## do lines intersect anywhere?
    denom=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    if denom*denom < epsilon:
        return False

    ## the lines do intersect, so let's see if they intersect
    ## inside both line segments
    t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4))/denom
    u = -1 * ((x1-x2)*(y1-y3) - (y1-y2)*(x1-x3))/denom

    ## do we care about falling inside the line segments? if so,
    ## check that the intersection falls within
    if inside and ( t < 0.0 or t > 1.0 or u < 0.0 or u > 1.0):
        return False

    return [x1 + t*(x2-x1), y1+t*(y2-y1), z1, 1.0]

## for a point and a line that lie in the x,y plane, compute the
## closest distance point on the line to the point, and return that
## point. If inside is true, then return the closest distance point
## between the point and the line segment.  If distance is true,
## return the distance, not the point.
##
## NOTE: in the case where inside is False, it's faster to compute the
## distance with linepoint(... inside=False, distance=True) than it is
## to compute the intersection point, so that is prefereable to the
## two-step operation of calling linePoint() and then dist() on the
## result.

def linePointXY(l,p,inside=True,distance=False):
    a=l[0]
    b=l[1]
    # check for degenerate case of zero-length line
    abdist = dist(a,b)
    if abdist < epsilon:
        return False

    x0=p[0]
    y0=p[1]
    z0=p[2]
    x1=a[0]
    y1=a[1]
    z1=a[2]
    x2=b[0]
    y2=b[1]
    z2=b[2]

    ## check to see if all three points lie in the same x,y plane
    if not isXYPlanar([p,a,b]):
        raise ValueError('non-XY points in linePointXY call')
        return false
    # if abs(z1-z0) > epsilon or abs(z2-z0) > epsilon:
    #    return False

    linedist = abs( ((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/abdist)

    ## this is the fast case:
    if not inside and distance:
        return linedist
    
    ## find out where the intersection between the original line and a
    ## line defined by the point and an orthogonal direction vector
    ## is.  We do this by constructing two direction vectors
    ## orthogonal to the orgiginal line scaled by the line distance,
    ## and adding them to the point in question.  Assuming that the
    ## line distance is not zero, only one of these constructed points
    ## will fall on the line

    ## compute unit direction vector for original line
    dir = sub(b,a)
    dir = scale(dir,1.0/mag(dir))

    ## compute two orthogonal direction vectors of length linedist
    ordir1 = scale(ortho(dir),linedist)
    ordir2 = scale(ordir1, -1.0)
    
    ## there are two possible intersection points
    pi1 = add(p,ordir1)
    pi2 = add(p,ordir2)

    ## compute distances
    d1 =  dist(a,pi1)+dist(pi1,b) # "triangle" with pi1
    d2 =  dist(a,pi2)+dist(pi2,b) # "triangle" with pi2

    ## the shortest "triangle" distance will signal the point that
    ## is actually on the line, even if that point falls outside
    ## the a,b line interval
    
    if not inside: # if we don't care about being inside the line
                   # segment
        if d1 <= d2:
            if distance:
                return d1
            else:
                return pi1
        else:
            if distance:
                return d2
            else:
                return pi2
        
    
    ## if the closest point on the line to point p lies between
    ## the endpoints of the line, then either d1 or d2 will equal
    ## abdist.  IF neither do, then we know that the closest point lies
    ## outside the endpoints

    if abs(d1-abdist) < epsilon:
        if distance:
            return linedist
        else:
            return pi1

    if  abs(d2-abdist) < epsilon:
        if distance:
            return linedist
        else:
            return pi2

    ## closest point is outside the interval.  That means that the
    ## distance from point p to whichever endpoint is smaller is the
    ## closest distance

    d3 = dist(a,p)
    d4 = dist(b,p)

    if d3 < d4:
        if distance:
            return d3
        else:
            return a
    else:
        if distance:
            return d4
        else:
            return b

## convenience function for fast distance calc
def linePointXYDist(l,p,inside=True):
    return linePointXY(l,p,inside,distance=True)
    
## functions that compute barycentric coordinates for 2D triangles and
## 3D tetrahedral simplexes

## given a point a and three vertices, all of which fall into the x-y
## plane, compute the barycentric coordinates lam1, lam2, and lam3
def barycentricXY(a,p1,p2,p3):
    if not isXYPlanar([a,p1,p2,p3]):
        raise ValueError('non-XY points in barycentricXY call')
    
    x=a[0]
    y=a[1]

    x1=p1[0]
    y1=p1[1]

    x2=p2[0]
    y2=p2[1]

    x3=p3[0]
    y3=p3[1]

    denom=(y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)

    if abs(denom) < epsilon:
        raise ValueError('degenerate triangle in barycentricXY call')
    
    lam1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/denom
    lam2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/denom
    lam3 = 1.0-(lam1+lam2)

    return [lam1,lam2,lam3]

## given a triangle defined as a list of three points and a test
## point, determine if that test point lies inside or outside the
## triangle, assuming that all points lie in the x-y plane

def _isinsidetriangleXY(a,poly): #no value checking version
    p1=poly[0]
    p2=poly[1]
    p3=poly[2]
    bary = barycentricXY(a,p1,p2,p3)
    if bary[0] >= 0.0 and bary[0] <= 1.0 and \
       bary[1] >= 0.0 and bary[1] <= 1.0 and \
       bary[2] >= 0.0 and bary[2] <= 1.0:
        return True
    return False
        

# 'safe' value-check version
def isinsidetriangleXY(a,poly):
    if len(poly) != 3 or not isXYPlanar(poly + [a ]):
        raise ValueError('bad poly length or non-XY points in insidetriangleXY call')
    return _isinsidetriangleXY(a,poly)

## given a polygon defined as a list of three or more points and a
## test point, determine if that test point lies inside or ouside the
## polygon, assuming that all points lie in the x-y plane

def _isinsidepolyXY(a,poly): #no value checking version
    if len(poly) == 3:
        return _isinsidetriangleXY(a,poly)
    else:
        return _isinsidetriangleXY(a,poly[0:3]) or \
               _isinsidepolyXY(a,poly[1:len(poly)])

# 'safe' value-check version
def isinsidepolyXY(a,poly):
    if len(poly) < 3 or not isXYPlanar(poly + [ a]):
        raise ValueError('bad poly length or non-XY points in insidepolyXY call')
    return _isinsidepolyXY(a,poly)

## check to see if we have been invoked on the command line
## if so, run some tests

if __name__ == "__main__":
    a = vect(5,0)
    b = vect(0,5)
    c = vect(-3,-3)
    d = vect(1,1)
    e = vect(10,10)
    
    print("light-weight unit tests for geom.py module")
    print("some points: a:" + vstr(a) + ", b:"
          + vstr(b) + ", c:" + vstr(c) + ", d:" + vstr(d) + ", e:" + vstr(e))

    print("mag a: " + str(mag(a)))
    print("add(a,b): " + vstr(add(a,b)))
    print("sub(a,b): " + vstr(sub(a,b)))
    print("mag(sub(a,b)): " + str(mag(sub(a,b))))
    print("mag(sub(a,b)) == sqrt(50): " + str(mag(sub(a,b))==sqrt(50.0)))
    
    l1 = [a,b]
    l2 = [c,d]
    l3 = [c,e]
    
    print("l1:" + vstr(l1) + ", l2:" + vstr(l2) +", l3:" + vstr(l3))

    int0 = intersectXY(l1,l1)
    int1 = intersectXY(l1,l2,False)
    int2 = intersectXY(l1,l2,True)
    int3 = intersectXY(l1,l3,True)

    print("expect False: intersectXY(l1,l1): " + vstr(int0))
    print("expect [2.5, 2.5]: intersectXY(l1,l2,False): " + vstr(int1))
    print("expect False: intersectXY(l1,l2,True): " + vstr(int2))
    print("expect [2.5, 2.5]: intersectXY(l1,l3,True): " + vstr(int3))
    
    print("linePointXY(l1,vect(0,0)): "
          + vstr(linePointXY(l1,vect(0,0))))
    print("linePointXYDist(l1,vect(0,0)) == sqrt(12.5): "
          + vstr(abs(linePointXYDist(l1,vect(0,0))-sqrt(12.5))<epsilon))
    print("linePointXY(l1,vect(0,10),False): "
          + vstr(linePointXY(l1,vect(vect(0,10)),False)))
    print("linePointXY(l1,vect(0,10),True): "
          + vstr(linePointXY(l1,vect(0,10),True)))
    print("linePointXY(l1,vect(10,0),False): "
          + vstr(linePointXY(l1,vect(10,0),False)))
    print("linePointXY(l1,vect(10,0),True): "
          + vstr(linePointXY(l1,vect(10,0),True)))

    print("do some planar point checking")
    # points in the x-y plane
    p1 = vect(2.5,0)
    p2 = vect(5,5)
    p3 = vect(0,5)
    p4 = vect(2.5,10)

    # points in the x,z plane
    p5 = vect(1,0,-5)
    p6 = vect(10,0,10)
    
    # points in the y-z plane
    p7 = vect(0,2.-1)
    p8 = vect(0,10,10)

    print('expect True: isCardinalPlanar("xy",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("xy",[p1,p2,p3,p4])))
    print('expect False: isCardinalPlanar("xz",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("xz",[p1,p2,p3,p4])))
    print('expect False: isCardinalPlanar("yz",[p1,p2,p3,p4]) : {}'.format(
        isCardinalPlanar("yz",[p1,p2,p3,p4])))

    print('expect True: isCardinalPlanar("xz",[p1,p5,p6]) : {}'.format(
        isCardinalPlanar("xz",[p1,p5,p6])))

    print('expect True: isCardinalPlanar("yz",[p2,p7,p8]) : {}'.format(
        isCardinalPlanar("yz",[p2,p7,p8])))

    try:
        print('deliberate bad plane specification, should raise ValueError')
        print('isCardinalPlanar("FOO!",[p1,p2,p3,p4]) : {}'.format(
            isCardinalPlanar("FOO!",[p1,p2,p3,p4])))
    except ValueError as err:
        print('got expected result:',err)

    print("polygon inside testing")
    tri1 = [p1,p2,p3]
    poly1 = [p1,p2,p4,p3]
    p = vect(2.5,1.0)
    q = vect(2.5,5.0)
    r = vect(2.5,7.0)
    s= vect(-10,-10)
    t= vect(10,10)
    print ("p: {}, q: {}, r: {}, s: {}, t: {}".format(vstr(p), vstr(q),
                                                      vstr(r), vstr(s),
                                                      vstr(t)))
    print ("tri1: {}".format(vstr(tri1)))
    print("expect True: isinsidetriangleXY(p,tri1): {}"\
          .format(isinsidetriangleXY(p,tri1)))
    print("expect True: isinsidetriangleXY(q,tri1): {}"\
          .format(isinsidetriangleXY(q,tri1)))
    print("expect False: isinsidetriangleXY(r,tri1): {}"\
          .format(isinsidetriangleXY(r,tri1)))
    print("expect False: isinsidetriangleXY(s,tri1): {}"\
          .format(isinsidetriangleXY(s,tri1)))
    print("expect False: isinsidetriangleXY(t,tri1): {}"\
          .format(isinsidetriangleXY(t,tri1)))
    print("inside poly testing")
    print ("tri1: {}".format(vstr(tri1)))
    print ("poly1: {}".format(vstr(poly1)))
    print("expect True: isinsidepolyXY(p,tri1): {}"\
          .format(isinsidepolyXY(p,tri1)))
    print("expect True: isinsidepolyXY(q,tri1): {}"\
          .format(isinsidepolyXY(q,tri1)))
    print("expect False: isinsidepolyXY(r,tri1): {}"\
          .format(isinsidepolyXY(r,tri1)))
    
    print("expect True: isinsidepolyXY(p,poly1): {}"\
          .format(isinsidepolyXY(p,poly1)))
    print("expect True: isinsidepolyXY(q,poly1): {}"\
          .format(isinsidepolyXY(q,poly1)))
    print("expect True: isinsidepolyXY(r,poly1): {}"\
          .format(isinsidepolyXY(r,poly1)))
    print("expect False: isinsidepolyXY(s,tri1): {}"\
          .format(isinsidepolyXY(s,tri1)))
    print("expect False: isinsidepolyXY(t,tri1): {}"\
          .format(isinsidepolyXY(t,tri1)))
    
    print("done!")

