# yapCAD
yet another procedural CAD and computational geometry system written in python 3

## goals

The purpose of yapCAD is to support 2D and 3D computational geometry and CAD projects in python3.  yapCAD is designed to support multiple rendering back-ends, such that a relatively small amount of code is necessary to add support for a 2D or 3D cad or drawing file format.

The foundations of yapCAD are grounded in decades of the author's experience with graphics system programming, 3D CAD and simulation. At the same time, yapCAD should make the easy stuff easy, and the more advanced stuff

The initial implementation of yapCAD provides DXF file creation support through the awesome ezdxf package.

## examples

It's pretty easy to make a DXF drawing with yapCAD.  Here is an example:

	from ezdxf_drawable import *
	from geom import *

	#set up DXF rendering
	drawable=ezdxfDraw()

    ## make some geometry

    # make a point located at 10,10 in the x-y plane, rendered as a small
    # cross and circle
    point=Point(vect(10,10),"xo")

    # make a line segment between the points -5,10 and 10,-5 in the x-y plane
    line=Line(vect(-5,10),
	          vect(10,-5))

    # make an arc with a center at 0,3 with a radius of 3, from 45 degrees
    # to 135 degrees
    arc=Arc(vect(0,3),3,45,135)

    # Render the geometry we've just made
    point.draw()
    line.draw()
    arc.draw()

    # write out the geometry as example1-out.dxf
    drawable.saveas("example1-out")

The yapCAD system isn't just about rendering, of course, it's about computational geometry.  For example, if you want to calculate the intersection of two lines in a plane, we have you covered

	from geom import *

    # define some points
    a = vect(5,0)
    b = vect(0,5)
    c = vect(-3,0)
    d = vect(10,10)

    # make a couple of lines
    l1 = [a,b]
    l2 = [c,d]

    # calculate the intersection of l1 and l2
	int0 = intersectXY(l1,l2)
	print("intersection of {} and {}: {}".format(l1,l2,int0))
	
## architecture

Under the hood, yapCAD is using [projective coordiates](https://en.wikipedia.org/wiki/Homogeneous_coordinates), sometimes called homogeneous coordinates, to represent points as 3D coodinates in the w=1 hyperplane. If that sounds complicated, its because it is. :P  But it does allow for a wide range of geometry operations, specifically [affine transforms](https://www.cs.utexas.edu/users/fussell/courses/cs384g-fall2011/lectures/lecture07-Affine.pdf) to be represented as composable transformation matricies. The benefits of this conceptual complexity is an architectual elegance and generality.

What does that buy you? It means that under the hood, yapCAD uses the same type of geometry engine that advanced CAD and GPU-based rendering systems use, and should allow for a wide range of computational geomety systems, possibly hardware-accelerated, to be built on top of it.

The good news is that you don't need to know about homogeneous coordinates, affine transforms, etc., to use yapCAD.  And most of the time you can pretend that your vectors are just two-dimensional if everything you are doing happens to lie in the x-y plane.

So, if you want to do simple 2D drawings, we have you covered.  If you want to buid a GPU-accelerated constructive solid geometry system, you can do that, too.
