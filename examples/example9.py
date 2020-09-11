## yapCAD 3D geometry example

print("example9.py -- yapCAD 3D geometry demonstration")
print("""

This is a demonstration of the creation and rendering of
three-dimensional geometry using yapCAD and the pyglet openGL draing
engine.

In this example we create an icosohedron centered at the orign and
tesellate it spherically""")

from geom import *
import math

## convert lat, lon, rad to cartesian coordinates
def lat_lon_rad(lat,lon,rad):
    if lat == 90:
        return [0,0,rad,1]
    elif lat == -90:
        return [0,0,-rad,1]
    else:
        latr = (((lat+90)%180)-90)*pi2/360.0
        lonr = (lon%360)*pi2/360.0
    
        smallrad = math.cos(latr)*rad
        z = math.sin(latr)*rad
        x = math.cos(lonr)*smallrad
        y = math.sin(lonr)*smallrad
        return [x,y,z,1]
    
def makeIcoPoints(center,radius):
    points = []
    normals = []
    p = lat_lon_rad(90,0,radius)
    n = scale(p,1.0/mag(p))
    points.append(p)
    normals.append(n)
    for i in range(10):
        sgn = 1
        if i%2 == 0:
            sgn =-1
        lat = math.atan(0.5)*360.0*sgn/pi2
        lon = i*36.0
        p = lat_lon_rad(lat,lon,radius)
        n = scale(p,1.0/mag(p))
        points.append(p)
        normals.append(n)

    p = lat_lon_rad(-90,0,radius)
    n = scale(p,1.0/mag(p))
    points.append(p)
    normals.append(n)
    return list(map( lambda x: add(x,center),points)),normals

icaIndices = [ [1,11,3],[3,11,5],[5,11,7],[7,11,9],[9,11,1],
               [2,1,3],[2,3,4],[4,3,5],[4,5,6],[6,5,7],[6,7,8],[8,7,9],[8,9,10],[10,9,1],[10,1,2],
               [0,2,4],[0,4,6],[0,6,8],[0,8,10],[0,10,2] ]


        
if __name__ == "__main__":
    verts,normals = makeIcoPoints(point(0,0,0),5.0)
    indices = icaIndices

    from pyglet_drawable import *
    dd=pygletDraw()
    dd.set_linecolor('white')

    for nd in indices:
        dd.draw(line(verts[nd[0]],
                     verts[nd[1]]))
        dd.draw(line(verts[nd[1]],
                     verts[nd[2]]))
        dd.draw(line(verts[nd[2]],
                     verts[nd[0]]))

    dd.draw_surface(verts,normals,indices)
    dd.set_linecolor('aqua')
    dd.polystyle = 'points'
    dd.pointstyle = 'xo'
    dd.draw(verts)
    dd.display()
