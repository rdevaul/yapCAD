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

def vect(a=False,b=False,c=False,d=False):
    r = [0,0,0,1]
    if isgoodnum(a):
        r[0]=a
        if isgoodnum(b):
            r[1]=b
            if isgoodnum(c):
                r[2]=c
                if isgoodnum(d):
                    r[3]=d
    elif isinstance(a,(tuple,list)):
        for i in range(min(4,len(a))):
            x=a[i]
            if isgoodnum(x):
                r[i]=x
    return r

## function to deep-copy geometry, which is to say lists containing
## non-zero-length lists or 4-vectors.  Returns False if argument
## isn't valid geometry list

def deepcopy(a):
    if isvect(a):
        return vect(a)
    elif isinstance(a,list):
        if len(a) > 0:
            c = list(map(deepcopy,a))
            for i in c:
                if not i:
                    return False
            return c
        else:
            return False
    else:
        return False

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
def orthoXY(a):
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
    
## operations on points
## --------------------

## points are defined as vectors that lie in a positive, non-zero
## hyperplane, i.e. [x, y, z, w] such that w > 0.  points are
## distinguished from pesuedovectors, such as the parameters to an
## arc, which don't transform as vectors.

## pseudovectors, by contrast, should lie in a w < 0 hyperplane.  For
## example, by convention arc parameters lie in the pseudovector w=-1
## hyperplane

## since points are vectors, we don't have a lot of special operations
## on points, except to explicitly create them and to test for them.

def point(x=False,y=False,z=False,w=False):
    if ispoint(x):
        return deepcopy(x)
    elif isgoodnum(x) and isgoodnum(y) and isgoodnum(z) and isgoodnum(w)\
       and w > 0:
        return [ x,y,z,w ]
    return vect(x,y,z)

def ispoint(x):
    if isvect(x) and x[3] > 0.0:
        return True
    return False
    
## operations on lines
## --------------------

## lines are defined as lists of two points, i.e.  [vect(x1,
## y1),vect(x2, y2)].  Lines must lie in a positive hyperplane,
## viz. w>0

## make a line, copying points, value-safe
def line(p1,p2=False):
    if isline(p1):
        return deepcopy(p1)
    elif ispoint(p1) and ispoint(p2):
        return [ point(p1), point(p2) ]
    else:
        raise ValueError('bad values passed to line()')

## is it a line?
def isline(l):
    return isinstance(l,list) and len(l) == 2 \
        and ispoint(l[0]) and ispoint(l[1])

## return the length of a line
def linelength(l):
    return dist(l[0],l[1])

## return the center of a line
def linecenter(l):
    return scale(add(l[0],l[1]),0.5)

## Sample a parameterized line.  Values 0 <= u <= 1.0 will fall within
## the line segment, values u < 0 and u > 1 will fall outside the line
## segment.

def sampleline(l,u):
    p1=l[0]
    p2=l[1]
    p = 1.0-u
    return add(scale(p1,p),scale(p2,u))

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
        raise ValueError('lines not in same x-y plane')

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
    ordir1 = scale(orthoXY(dir),linedist)
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


## functions that operate on 2D arcs/circles

## an arc or circle is defined by a center, a radius, a start angle,
## an end angle, and a normal that specifies the plane of the
## arc/cicle.

## angles are in degress, and specify a counter-clockwise sweep.  If
## start = 0 and end = 360 (both integer) the arc is a circle.  Other
## values that create a 360 degree difference may produce an
## infentesimal gap.

## an arc is defined as a list of one or two vectors and a
## pseudovector: [ center, [r, s, e], <normal>].  Most of the time, we
## assume that arcs lie in the x-y plane, which is to say that
## normal = [0,0,1] if it's not specified.

## to mark that the second list element is a pseudovector, we set the
## w component to -1.  Since negative w values should never exist for
## vectors in our projective geometry system this should be a robust
## convention.

## make an arc, copying points, value-safe
## NOTE: if start and end are not specified, a full circle is created
def arc(c,rp=False,sn=False,e=False,n=False):
    if isarc(c):
        return deepcopy(c)
    elif ispoint(c):
        cen = point(c)
        if isvect(rp):
            psu = deepcopy(rp)
            r=psu[0]
            if r < 0:
                raise ValueError('negative radius not allowed for arc')
            psu[3]=-1
            if not sn:
                return [ cen, psu ]
            else:
                if not ispoint(sn) or abs(mag(sn)-1.0) > epsilon:
                    raise ValueError('bad (non-unitary) plane vector for arc')
                return [ cen, psu, point(sn) ]
        elif isgoodnum(rp):
            r = rp
            if r < 0:
                raise ValueError('negative radius not allowed for arc')
            start = 0
            end = 360
            if isgoodnum(sn) and isgoodnum(e):
                start = sn
                end = e
            psu = vect(r,start,end,-1)
            if not n:
                return [ cen, psu ]
            else:
                if not ispoint(n) or abs(mag(n)-1.0) > epsilon:
                    raise ValueError('bad (non-unitary) plane vector for arc')
                return [ cen, psu, point(n) ]

    raise ValueError('bad arguments passed to arc()')
        
## function to determine if argument can be interpreted as a valid
## arc.

def isarc(a):
    if not isinstance(a,list):
        return False
    n = len(a)
    if n < 2 or n > 3:
        return False
    if not (ispoint(a[0]) and isvect(a[1])):
        return False
    if a[1][3] != -1:           # is psuedovector marked?
        return False
    r =a[1][0]                
    if r < 0:
        # is radius less than zero? if so, not a valid arc
        return False
    if n == 3 and ( not ispoint(a[2]) or abs(mag(a[2])-1.0) > epsilon):
        # if plane-definition vector is included but is non-unitary,
        # it's not a valid arc
        return False
    
    return True

## test to see if an arc is a circle.  NOTE that due to modulus
## arithmetic, it's not possible to specify an arc with a full 360
## degees range, as this would "wrap around" to an arc of 0 degrees
## range.  To solve this, we use a special convention for start and
## end to signal a true full circle

def iscircle(a):
    if isarc(a):
        start=a[1][1] 
        end=a[1][2]
        ## these are special, integer values that flag a true full
        ## circle.
        if start==0 and end==360:
            return True
    else:
        return False

## function to return the midpoint of an arc
def arccenter(c):
    return samplearc(c,0.5)
    
## function to return the length of an arc
def arclength(c):
    r=c[1][0]
    start=c[1][1] % 360.0
    end=c[1][2]%360.0

    d = pi2*r
    l = d*(start-end)/360.0
    return l

def samplearc(c,u):
    p=c[0]
    r=c[1][0]
    start=c[1][1] % 360.0
    end=c[1][2] %360.0

    if len(c) == 3:
        norm = c[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise NotImplementedError('non x-y plane arc sampling not yet supported')
    angle = ((end-start)*u+start)
    radians = angle*pi2/360.0
    q = scale(vect(cos(radians),sin(radians)),r)

    return add(p,q)

## Intersection functions for arcs and circles.
## we can compute the intersection of an arc and a line, or an arc and
## an arc, when these all lie in the same plane.

## arc-arc jtersection calculation, non-value-safe version
def _arcArcIntersectXY(c1,c2,inside=True):
    x1=c1[0]
    x2=c2[0]
    r1=c1[1][0]
    r2=c2[1][0]

    ## first check for non-intersection due to distance between the
    ## centers of the arcs, treating both arcs as circles for the moment

    d=dist(x1,x2) #calculate the distance d between circle centers

    if d > r1+r2:
        return False # too far, no possiblity of intersection

    if ( r1> r2 and d < r1-r2) or (r2 >= r1 and d < r2-r1):
        return False # too close, little arc is fully inside bigger arc

    if d < epsilon:
        return False # circle centers too close for stable calculation

    ## OK, we are in the goldilocks zone of intersection.  this means
    ## that if boh arcs are cicles or if inside=False we are
    ## guaranteed one or two intersections.  Calculate those
    ## intersections and then test to see if they fall between start
    ## and end of the respective arcs

    ## we start by calculating the distance id of the intersection plane
    ## from the center of arc 1, knowing that by definition id <= r1

    ## Math: consider the triangle with side lengths r1, r2, and d,
    ## where d is the previously calculated distance between arc
    ## centers.  Consider the two right triangles with side lengths
    ## r1, id, h, and r2, h, (d-id).  We know that:
    ## id^2 + h^2 = r1^2, (d-id)^2 + h^2 = r2^2
    ## solving both for h2 and then substituting, this means:
    ## r1^2 - id^2 = r2^2 - (d-id)^2
    ## collecting terms and solving for id produces:
    ## id = (r1^2-r2^2 + d^2)/2d

    id = (r1*r1 - r2*r2 + d*d)/(2 * d)

    ## compute the point on the line connecting the two arc centers
    ## that is id away from the first arc

    v1 = scale(sub(x2,x1),1.0/d) # unitary direction vector pointing
                                 # from x1 to x2
    v2 = scale(v1,id) # point on line between two circles in
                      # coordinate space centered at x1

    ## compute direction vector o orthgonal to v1 -- the line that
    ## intersects point v2 and v2+o will pass through our intersection
    ## points

    o = orthoXY(v1)
    
    ## now, transform back into world coordinates and calculate the
    ## intersection of this line with either of our arcs, treating
    ## them as circles for now

    l = [add(v2,x1),add(add(v2,o),x1)]

    s = _lineArcIntersectXY(l,c1,False)

    ## as a sanity check, do the same with the other arc.  Results
    ## should be within epsilon
    #ss = _lineArcIntersectXY(l,c2,False)
    #foo = list(map(lambda x, y: dist(x,y) < epsilon,s,ss))
    #print("sanity check: " , foo)

    if not s or len(s) == 0:
        raise ValueError('no computed intersections, something is wrong')

    if not inside:
        return s
    
    ## jump back to arc1 and arc2 space and check angles

    s1 = list(map(lambda x: sub(x,x1),s))
    s2 = list(map(lambda x: sub(x,x2),s))

    ## compute start and end angles for arcs
    start1=c1[1][1] % 360.0
    end1=c1[1][2] %360.0

    start2=c2[1][1] % 360.0
    end2=c2[1][2] %360.0

    ## check each intersection against angles for each arc.  
    ss = []
    for i in range(len(s)):
        p1 =s1[i]
        p2 =s2[i]
        ang1 = (atan2(p1[1],p1[0]) % pi2)*360.0/pi2
        ang2 = (atan2(p2[1],p2[0]) % pi2)*360.0/pi2

        good = False
        ## check angle against first arc
        if end1 > start1 and ang1 >= start1 and ang1 <= end1:
            good = True
        elif start1 < end1 and ang1 <= start1 and ang1>= end1:
            good = True

        ## check angle against second arc
        if end2 > start2 and ang2 >= start2 and ang2 <= end2:
            good = good and True
        elif start2 < end2 and ang2 <= start2 and ang2>= end2:
            good = good and True
        else:
            good = False

        ## only add instersection to the list if both checks were passed
        if good:
            ss = ss + [ s[i] ]
    if len(ss) == 0:
        return False
    else:
        return ss

## value-safe wrapper for arc-arc intersection function
def arcArcIntersectXY(c1,c2,inside=True):
    for c in [c1,c2]:
        if len(c) == 3:
            norm = c[2]
            if dist(norm,vect(0,0,1)) > epsilon:
                raise ValueError('arc passed to lineArcIntersectXY does not lie in x-y plane')
    if not isXYPlanar([c1[0],c2[0]]):
        raise ValueError('arcs passed to arcArcIntersectXY do not lie in same x-y plane')
    return _arcArcIntersectXY(c1,c2,inside)
    

## non-value-safe line-arc intersection function
def _lineArcIntersectXY(l,c,inside=True):
    x=c[0]
    r=c[1][0]

    # is the arc a full circle?
    circle = False
    if c[1][1] == 0 and c[1][2] == 360:
        circle = True
        
    start=c[1][1] % 360.0
    end=c[1][2] %360.0

    ## what is the shortest distance between the line and the center
    ## of the arc?  If that is greater than r, then there is no
    ## intersection
    dist = linePointXYDist(l,x,inside)
    if dist > r:
        return False

    ## start by treating the arc as a circle.  At this point we know
    ## we have one or two intersections within the line segment,
    ## though perhaps none within the arc segment, which we will test
    ## for later
    
    ## transform points so arc is located at the origin
    p0=sub(l[0],x)
    p1=sub(l[1],x)
    
    ## solve for b in:  | b*p0 + (1-b)*p1 | = r
    ## let V= p0-p1, P=p1
    ##     | b*V + P |^2 = r^2
    ##       b^2(Vx^2 + Vy^2) + 2b(VxPx+VyPy) + Px^2 + Py^2 - r^2 = 0
    ## let a = Vx^2 + Vy^2,
    ##     b = 2*(VxPx + VyPy)
    ##     c = Px^2 + Py^2 - r^2
    ## b0 = ( -b + sqrt(b^2 - 4ac) )/ 2a
    ## b1 = ( -b - sqrt(b^2 - 4ac) )/ 2a
    
    V = sub(p0,p1)
    P = p1
    a = V[0]*V[0]+V[1]*V[1]
    if abs(a) < epsilon:
        raise ValueError('degenerate line in lineArcIntersectXY')
    b = 2*(V[0]*P[0]+V[1]*P[1])
    c = P[0]*P[0]+P[1]*P[1]-r*r
    d = b*b-4*a*c
    if d < 0:
        raise ValueError("imaginary solution to circle line intersection -- shouldn't happen here")
    if d < epsilon: # one point of intersection
        b0 = -b/(2*a)
        b1 = False
    else: # two points of intersection
        b0 = (-b + sqrt(d))/(2*a)
        b1 = (-b - sqrt(d))/(2*a)

    # use computed parameters to calculate solutions, still in
    # circle-at-origin coordinates
    s = [ add(scale(V,b0),p1) ]
    if b1:
        s = s + [ add(scale(V,b1),p1) ]

    if not inside or circle:              # transform back into world
                                          # coordinates
        return list(map(lambda q: add(q,x),s))

    ## see if any of the intersections we've found lie between
    ## start and end of the arc
    
    ss = []
    for i in s:
        ang = (atan2(i[1],i[0]) % pi2)*360.0/pi2

        if end > start and ang >= start and ang <= end:
            ss = ss + [ add(x,i) ]
        elif start < end and ang <= start and ang>= end:
            ss = ss + [ add(x,i) ]

    if len(ss) == 0:
        return False
    return ss

## value-safe wrapper for line-arc intersection function
def lineArcIntersectXY(l,c,inside=True):
    if len(c) == 3:
        norm = c[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('arc passed to lineArcIntersectXY does not lie in x-y plane')
    points = l + [ c[0] ]
    if not isXYPlanar(points):
        raise ValueError('line and circle passed to lineArcIntersectXY do not all lie in same x-y plane')
    return _lineArcIntersectXY(l,c,inside)


## function to compute tangent lines to two coplanar circles lying in
## an x-y plane.  Function will either return two lines or False, if
## the center of the circles are too close together.

def _circleCircleTangentsXY(c1,c2):
    a = c1[1][0]
    b = c2[1][0]
    if a>b:
        bigIsOne=True
        bigC = c1
        smallC = c2
    else:
        bigIsOne=False
        bigC = c2
        smallC = c1
    ## Consdier the triangle created by the center of the small
    ## circle, the center of the large circle, and the point at the 90
    ## degree intersection of the line from the center of the small
    ## circle to the radian of the tangent point on the large circle.
    ## This is a right triangle with one leg of length d (distance of
    ## centers), one leg of length bigR-smallR, and one leg of unknown
    ## length, beta. theta is the angle formed by d and beta, which is
    ## also the angle of one of the the tangent lines, the other being
    ## -theta.
    ## 
    ## we will calulate theta as follows:
    ## beta^2 - (r2-r1)^2 = d^2
    ## beta = sqrt( d^2 - (r2-r1)^2 )
    ## theta = atan ((r2-r1)/beta)
    
    r1 = smallC[1][0]
    r2 = bigC[1][0]

    d = dist(c1[0],c2[0])
    dr = r2-r1

    if d <= dr: #centers too close
        return False
    
    beta = sqrt( d*d - dr*dr)
    theta = atan2(dr,beta)

    ## now, figure out the angle created by the center of the large
    ## circle with respect to the small circle
    dd = sub(bigC[0],smallC[0])
    phi = atan2(dd[1],dd[0])

    ## the two lines have angle phi+theta, and phi-theta.  The
    ## intersection point of these lines is at the point on the circle
    ## phi+theta+90', and phi-theta-90'
    gamma1 = phi+theta+pi/2
    gamma2 = phi-theta-pi/2
    n1 = point(cos(gamma1),sin(gamma1))
    n2 = point(cos(gamma2),sin(gamma2))
    p1 = add(scale(n1,r1),smallC[0])
    p2 = add(scale(n1,r2),bigC[0])
    p3 = add(scale(n2,r1),smallC[0])
    p4 = add(scale(n2,r2),bigC[0])

    l1 = line(p1,p2)
    l2 = line(p3,p4)

    return [l1,l2]

def circleCircleTangentsXY(c0,c1):
    if not iscircle(c0) or not iscircle(c1):
        raise ValueError('non circles passed to circleCirlceTangentsXY')
    if len(c0) == 3:
        norm = c0[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('first circle passed to circleCircleTangentsXY does not lie in an x-y plane')
    if len(c1) == 3:
        norm = c1[2]
        if dist(norm,vect(0,0,1)) > epsilon:
            raise ValueError('second circle passed to circleCircleTangentsXY does not lie in an x-y plane')
    points = [ c0[0], c1[0]]
    if not isXYPlanar(points):
        raise ValueError('circles passed to circleCircleTangentsXY do not all lie in same x-y plane')
    return _circleCircleTangentsXY(c0,c1)


## functions that compute barycentric coordinates for 2D triangles

## given a point a and three vertices, all of which fall into the x-y
## plane, compute the barycentric coordinates lam1, lam2, and lam3
def _barycentricXY(a,p1,p2,p3):
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

def barycentricXY(a,p1,p2,p3):
    if not isXYPlanar([a,p1,p2,p3]):
        raise ValueError('non-XY points in barycentricXY call')
    return _barycentricXY(a,p1,p2,p3)

## given a triangle defined as a list of three points and an
## additional test point, determine if that test point lies inside or
## outside the triangle, assuming that all points lie in the x-y plane

def _isInsideTriangleXY(a,poly): #no value checking version
    p1=poly[0]
    p2=poly[1]
    p3=poly[2]
    bary = _barycentricXY(a,p1,p2,p3)
    if bary[0] >= 0.0 and bary[0] <= 1.0 and \
       bary[1] >= 0.0 and bary[1] <= 1.0 and \
       bary[2] >= 0.0 and bary[2] <= 1.0:
        return True
    return False
        

# value-safe wrapper for triangle inside testing
def isInsideTriangleXY(a,poly):
    if len(poly) != 3 or not isXYPlanar(poly + [a ]):
        raise ValueError('bad poly length or non-XY points in insidetriangleXY call')
    return _isInsideTriangleXY(a,poly)

## given a convex polygon defined as a list of three or more points
## and a test point, determine if that test point lies inside or
## ouside the polygon, assuming that all points lie in the x-y plane

def _isInsideConvexPolyXY(a,poly): #no value checking version
    if len(poly) == 3:
        return _isInsideTriangleXY(a,poly)
    else:
        return _isInsideTriangleXY(a,poly[0:3]) or \
               _isInsideConvexPolyXY(a,poly[1:len(poly)])

# value-safe wrapper for convex polygon inside testing function
def isInsideConvexPolyXY(a,poly):
    if len(poly) < 3 or not isXYPlanar(poly + [ a]):
        raise ValueError('bad poly length or non-XY points in insidepolyXY call')
    return _isInsideConvexPolyXY(a,poly)


# -------------------------------------------------------------
## check to see if we have been invoked on the command line
## if so, run some tests

if __name__ == "__main__":
    print("------------------------------------------")
    print("yapCAD geometry tests")
    print("------------------------------------------")
    print("light-weight unit tests for geom.py module")
    a = point(5,0)
    b = point(0,5)
    c = point(-3,-3)
    d = point(1,1)
    e = point(10,10)
    print("---> point creation and testing")
    print("some points: a:" + vstr(a) + ", b:"
          + vstr(b) + ", c:" + vstr(c) + ", d:" + vstr(d) + ", e:" + vstr(e))
    print("ispoint(a): ",ispoint(a))
    print("ispoint(point(a)): ",ispoint(point(a)))
    print("ispoint([1,2]: ",ispoint([1,2]))
    print("ispoint(vect(1,2)): ",ispoint(vect(1,2)))
    print("ispoint(vect(1,2,3,4)): ",ispoint(vect(1,2,3,4)))
    print("ispoint(vect(1,2,3,-1)): ",ispoint(vect(1,2,3,-1)))
    
    print("---> basic vector operations tests")
    print("mag a: " + str(mag(a)))
    print("add(a,b): " + vstr(add(a,b)))
    print("sub(a,b): " + vstr(sub(a,b)))
    print("mag(sub(a,b)): " + str(mag(sub(a,b))))
    print("mag(sub(a,b)) == sqrt(50): " + str(mag(sub(a,b))==sqrt(50.0)))

    print("---> line creation and testing")
    l1 = [a,b]
    print("l1 = [a,b] -- l1:",vstr(l1))
    l11 = line(a,b)
    print("l11 = line(a,b) -- l1:",vstr(l11))
    l2 = [c,d]
    l22 = line(l2)
    print("l2 = [c,d], l22 = line(l2) -- l22: ",vstr(l22))
    l3 = line(c,e)
    print("l3 = line(c,e), isline(l3) : ",isline(l3))
    print("a {}, isline(a): {}".format(vstr(a),isline(a)))
    
    
    print("---> vector and geometry copying tests")
    foo = [a,b,l3,l2,d]
    print("foo: ",vstr(foo))
    print("deepcopy(foo)",vstr(deepcopy(foo)))
    bar = [a,b,[1,2],l2,l3]
    print("bar: ",vstr(bar))
    print("expect False: deepcopy(bar)",vstr(deepcopy(bar)))
        
    print("---> line-line intersection tests")
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

    print("---> arc creation and testing")
    arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0,-1)]
    print("arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0,-1)], arc1: ",vstr(arc1))
    arc11=arc(vect(2.5,2.5),2.5,90.0,270.0)
    print("arc11=arc(vect(2.5,2.5),2.5,90.0,270.0), arc11: ",vstr(arc11))
    print("isarc(arc1): {}  isarc(arc11): {}".format(isarc(arc1),isarc(arc11)))
    arc12=arc(arc11)
    print("arc12=arc(arc11), arc12: {}, isarc(arc12): {}".format(vstr(arc12),isarc(arc12)))
    try:
        print("try creating an arc with a negative radius, should raise ValueError")
        print("arc(vect(0,0),-2): ",arc(vect(0,0),-2))
    except ValueError as err:
        print('got expected result:',err)
    print("--> line-arc disambiguation")
    print("l1: ",vstr(l1)," arc1: ",vstr(arc1))
    print("isline(l1): {} isline(arc1): {}".format(isline(l1),isline(arc1)))
    print("isarc(l1): {} isarc(arc1): {}".format(isarc(l1),isarc(arc1)))
    print("---> arc-line intersection tests")
    arc1=[vect(2.5,2.5),vect(2.5,90.0,270.0)]
    print("arc1: {}".format(vstr(arc1)))
    print("l1: {}".format(vstr(l1)))
    l2[1]=vect(0,0)
    print("l2: {}".format(vstr(l2)))
    int4 = lineArcIntersectXY(l1,arc1,False)
    int5 = lineArcIntersectXY(l1,arc1,True)
    int6 = lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True)
    int7 = lineArcIntersectXY(l2,arc1,True)
    int8 = lineArcIntersectXY(l2,arc1,False)
    print("lineArcIntersectXY(l1,arc1,False): {}".format(vstr(int4)))
    print("lineArcIntersectXY(l1,arc1,True): {}".format(vstr(int5)))
    print("lineArcIntersectXY([vect(0,5),vect(5,5)],arc1,True): {}".format(vstr(int6)))
    print("lineArcIntersectXY(l2,arc1,False): {}".format(vstr(int7)))
    print("lineArcIntersectXY(l2,arc1,True): {}".format(vstr(int8)))

    print("---> circle-circle tangent testing")
    circ1 = arc(point(5,5),5)
    circ2 = arc(point(-5,5),7.5)
    circ3 = arc(point(0,0),1)

    tl1 = circleCircleTangentsXY(circ1,circ2)
    tl2 = circleCircleTangentsXY(circ2,circ1)
    tl3 = circleCircleTangentsXY(circ3,circ2)

    print("circ1: ",vstr(circ1))
    print("circ2: ",vstr(circ2))
    print("circ3: ",vstr(circ3))
    
    print("circleCircleTangentsXY(circ1,circ2) :", vstr(tl1))
    print("circleCircleTangentsXY(circ2,circ1) :", vstr(tl2))
    print("circleCircleTangentsXY(circ3,circ2) :", vstr(tl3))

    print("---> arc-arc intersection tests")
    arc2=[vect(4.0,2.5),vect(2.5,90.0,270.0)]
    print("arc1: {}".format(vstr(arc1)))
    print("arc2: {}".format(vstr(arc2)))

    int9 = arcArcIntersectXY(arc1,arc2,False)
    int10 = arcArcIntersectXY(arc1,arc2,True)
    int11 = arcArcIntersectXY(arc2,arc1,True)
    print("arcArcIntersectXY(arc1,arc2,False):",vstr(int9))
    print("arcArcIntersectXY(arc1,arc2,True):",vstr(int10))
    print("arcArcIntersectXY(arc2,arc1,True):",vstr(int11))

    
    print("---> do some planar point testing")
    # points in the x-y plane
    p1 = point(2.5,0)
    p2 = point(5,5)
    p3 = point(0,5)
    p4 = point(2.5,10)

    # points in the x,z plane
    p5 = point(1,0,-5)
    p6 = point(10,0,10)
    
    # points in the y-z plane
    p7 = point(0,2.-1)
    p8 = point(0,10,10)

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


    print("---> convex polygon inside testing")
    tri1 = [p1,p2,p3]
    poly1 = [p1,p2,p4,p3]
    p = point(2.5,1.0)
    q = point(2.5,5.0)
    r = point(2.5,7.0)
    s= point(-10,-10)
    t= point(10,10)
    print ("p: {}, q: {}, r: {}, s: {}, t: {}".format(vstr(p), vstr(q),
                                                      vstr(r), vstr(s),
                                                      vstr(t)))
    print ("tri1: {}".format(vstr(tri1)))
    print("expect True: isInsideTriangleXY(p,tri1): {}"\
          .format(isInsideTriangleXY(p,tri1)))
    print("expect True: isInsideTriangleXY(q,tri1): {}"\
          .format(isInsideTriangleXY(q,tri1)))
    print("expect False: isInsideTriangleXY(r,tri1): {}"\
          .format(isInsideTriangleXY(r,tri1)))
    print("expect False: isInsideTriangleXY(s,tri1): {}"\
          .format(isInsideTriangleXY(s,tri1)))
    print("expect False: isInsideTriangleXY(t,tri1): {}"\
          .format(isInsideTriangleXY(t,tri1)))
    print("inside poly testing")
    print ("tri1: {}".format(vstr(tri1)))
    print ("poly1: {}".format(vstr(poly1)))
    print("expect True: isInsideConvexPolyXY(p,tri1): {}"\
          .format(isInsideConvexPolyXY(p,tri1)))
    print("expect True: isInsideConvexPolyXY(q,tri1): {}"\
          .format(isInsideConvexPolyXY(q,tri1)))
    print("expect False: isInsideConvexPolyXY(r,tri1): {}"\
          .format(isInsideConvexPolyXY(r,tri1)))
    
    print("expect True: isInsideConvexPolyXY(p,poly1): {}"\
          .format(isInsideConvexPolyXY(p,poly1)))
    print("expect True: isInsideConvexPolyXY(q,poly1): {}"\
          .format(isInsideConvexPolyXY(q,poly1)))
    print("expect True: isInsideConvexPolyXY(r,poly1): {}"\
          .format(isInsideConvexPolyXY(r,poly1)))
    print("expect False: isInsideConvexPolyXY(s,tri1): {}"\
          .format(isInsideConvexPolyXY(s,tri1)))
    print("expect False: isInsideConvexPolyXY(t,tri1): {}"\
          .format(isInsideConvexPolyXY(t,tri1)))
    
    print("done!")

