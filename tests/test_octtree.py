import pytest
from yapcad.geom import *
from yapcad.octtree import *
import examples.example10 as example10

class TestBuild:
    """ Test octtree and quadtree construction """

    p0 = point(-5,-5,-5)
    p1 = point(5,5,5)
    l = line(p0,p1)

    rp = example10.randomPoints(l,20)
    rp2 = example10.randomPoints(l,20)
    a = arc(point(0,0,0),2)
    a2 = arc(point(10,0,0),2)

    gl = [ p0, p1, l, a, a2]
    
    foo =NTree()
    foo.addElement(l)

    for i in range(len(rp)):
        foo.addElement(line(rp[i],rp2[i]))
    
    foo.addElement(a)
    foo.addElement(a2)

    foo.updateTree()
    assert False
    box = bbox(gl)

    # foo.getElements(bbox)
