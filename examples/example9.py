## yapCAD 3D geometry example

from yapcad.geom import *
from yapcad.geom3d import *
from yapcad.geom3d_util import *

import math


if __name__ == "__main__":
    from yapcad.pyglet_drawable import *
    print("example9.py -- yapCAD 3D geometry demonstration")
    print("""
This is a demonstration of the creation and rendering of
three-dimensional geometry using yapCAD and the pyglet openGL draing
engine.

In this example we create an icosohedron centered at the orign and
tesellate it spherically""")

    dd=pygletDraw()

    # rotate the camera down
    dd.ry = -80.0

    # make a spherical surface
    sphereobj = sphere(50.0,point(0,0,0),3)

    # make a rectangular prism
    pris = prism(50,30,20)

    # make a cone

    cone = conic(20,4,20)

    # import pdb ; pdb.set_trace()
    # make a jade-like material
    materials['jade'] = deepcopy(materials['default'])
    materials['jade'].ambient = [0.2, 0.7, 0.4, 0.8]
    materials['jade'].diffuse = [0.3, 0.5, 0.3, 0.8]
    materials['jade'].shininess= 128

    
    # set up some parameters for animation

    # these timers could be consolidated
    time1 = 0
    time2 = 0

    # paths for objects to move on
    circle1 = rotate(geomlist2poly([ arc(point(-30,0,0),25) ]),
                     45,axis=point(1.0/sqrt(2),1.0/sqrt(2),0))
    circle2 = arc(point(30,0,0),50)
    
    def animateSphere1(dt):
        global time1
        time1 += dt
        obj = dd.objectdict['sphere1']

        pos = sample(circle1,(time1/10)%1.0)
        obj.x = pos[0]
        obj.y = pos[1]
        obj.z = pos[2]
        obj.rz = time1*10.0

    def animateSphere2(dt):
        global time2
        time2 += dt
        obj = dd.objectdict['sphere2']

        pos = sample(circle2,(time2/10+0.5)%1.0)
        obj.x = pos[0]
        obj.y = pos[1]
        obj.z = pos[2]
        obj.rz = -time2*10.0

    ## create some aminatable, renderable objects

    dd.make_object('sphere0',lighting=False,linecolor='white',
                   position=point(0,0,0))
    
    dd.make_object('sphere1',lighting=True,material='jade',
                   linecolor='green',
                   animate=animateSphere1)
    
    dd.make_object('sphere2',lighting=True,material='default',
                   linecolor='yellow',
                   animate=animateSphere2)

    dd.draw(sphereobj,name='sphere0')
    dd.draw(sphereobj,name='sphere1')
    #dd.draw(sphereobj,name='sphere2')
    #dd.draw(pris,name='sphere2')
    dd.draw(cone,name='sphere2')


    # for the stationary sphere, draw some additional points
    dd.linecolor = 'aqua'
    dd.polystyle = 'points'
    dd.pointstyle = 'xo'
    
    spheresurf = sphereobj[1][0]
    dd.draw(spheresurf[1]) # the verticies of the sphere
    dd.display()
