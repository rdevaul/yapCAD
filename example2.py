## First computational geometry example for yapCAD

from geom import *

# define some points
a = vect(5,0)
b = vect(0,5)
c = vect(-3,-3)
d = vect(10,10)

# make a couple of lines
l1 = [a,b]
l2 = [c,d]

# calculate the intersection of l1 and l2
int0 = intersectXY(l1,l2,True)
print("intersection of {} and {}: {}".format(l1,l2,int0))
