import pytest
from yapcad.geom import *
from yapcad.geom_util import *
from yapcad.geom3d import *
from yapcad.pyglet_drawable import *
from yapcad.poly import *
from yapcad.combine import *

class TestSurface:
    """test surface representation"""
    def test_surface(self):
        
        dd = pygletDraw()
        dd.linecolor='silver'
        # make an arc-segment geometry list spiral
        gl = makeArcSpiral(point(-10,0),10.0,1.0)
        # "close" the arc with a line segment
        gl.append(line(sample(gl,1.0),
                       sample(gl,0.0)))
        dd.draw(gl)
        # turn the geometry list into a poly

        ply = geomlist2poly(gl)
        ply = translate(ply,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply)
        surf,bnd = poly2surface(ply)
        dd.draw_surface(surf)

        ply2 = [ point(20,20),
                 point(10,20),
                 point(10,-20),
                 point(20,-20),
                 point(20,-10),
                 point(15,-10),
                 point(15,10),
                 point(20,10),
                 point(20,20) ]
                 
        ply2 = translate(ply2,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply2)
        surf2,bnd = poly2surface(ply2)
        dd.draw_surface(surf2)

        rr1 = makeRoundRect(5,5,0.5,center=point(-15,20))
        rr2 = makeRoundRect(6,3,0.5,center=point(-12,20))
        gl3 = Boolean('difference',[rr1,rr2])
        ply3 = geomlist2poly(gl3.geom())
        dd.linecolor='aqua'
        dd.draw(gl3)
        ply3 = translate(ply3,point(0,0,0.5)) # pull it forward
        surf3,bnd = poly2surface(ply3)
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply3)
        dd.draw_surface(surf3)

        gl4 = [arc(point(0,7),10.0,0.0,270.0),
               arc(point(0,7),5.0,0.0,270.0,samplereverse=True)]
        gl4 = [arc(point(0,7),10.0,0.0,350.0),
               arc(point(0,7),5.0,0.0,350.0,samplereverse=True)]
        dd.linecolor='aqua'
        dd.draw(gl4)
        ply4 =geomlist2poly(gl4,minang = 10 )
        ply4.append(ply4[0])

        ply4 = translate(ply4,point(0,0,0.5)) # pull it forward
        dd.linecolor='red'
        dd.polystyle='both'
        dd.draw(ply4)

        surf4,bnd = poly2surface(ply4,minlen=1.0)
        # assert False
        lines = surf2lines(surf4)
        lines = translate(lines,point(0,0,0.5))
        dd.draw_surface(surf4)
        dd.linecolor='yellow'
        dd.draw(lines)
        bnd.append(bnd[0])
        bnd = translate(bnd,point(0,0,1.0))
        dd.linecolor='white'
        dd.draw(bnd)
        print(f'len(surf4): {len(surf4)}')
        print(f'len vertices: {len(surf4[1])}')
        print(f'len normals: {len(surf4[2])}')
        print(f'len faces: {len(surf4[3])}')
        print(f'len boundaries: {len(surf4[4])}')
        print(f'len holes: {len(surf4[5])}')
        assert issurface(surf4)
        assert issurface(surf4,fast=False)
        dd.display()
