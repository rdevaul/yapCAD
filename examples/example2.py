## First computational geometry example for yapCAD
print("example2.py -- yapCAD computational geometry example")
print('''
In this example, we create points, lines and arcs, and compute the
intersection of lines and arcs.''')

from yapcad.geom import *

# define some points
a = point(5,0)
b = point(0,5)
c = point(-3,-3)
d = point(10,10)

# make a couple of lines
l1 = line(a,b)
l2 = line(c,d)

# calculate the intersection of l1 and l2
int0 = intersectXY(l1,l2,True)

# define a semicircular arc centerd at 2.5, 2,5 with a radius of 2.5
# extending from 90 degress to 135 degrees

arc1=arc(point(2.5,2.5),2.5,90.0,270.0)

# calculate the intersection of the line l1 and the arc arc1

int1 = intersectXY(l1,arc1,True)

# print the results
print("line l1:",vstr(l1))
print("line l2:",vstr(l2))      
print("arc arc1:",vstr(arc1))

print("intersection of l1 and l2:", vstr(int0))
print("intersection of arc1 and l1:",vstr(int1))

