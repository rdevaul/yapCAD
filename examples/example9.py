## yapCAD 3D geometry example

from yapcad.geom import *
import math

## utility function, convert lat, lon, rad to cartesian coordinates
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

## icosohedron-specific function, generate intial geometry
def makeIcoPoints(center,radius):
    points = []
    normals = []
    p = lat_lon_rad(90,0,radius)
    n = scale3(p,1.0/mag(p))
    points.append(p)
    normals.append(n)
    for i in range(10):
        sgn = 1
        if i%2 == 0:
            sgn =-1
        lat = math.atan(0.5)*360.0*sgn/pi2
        lon = i*36.0
        p = lat_lon_rad(lat,lon,radius)
        n = scale3(p,1.0/mag(p))
        points.append(p)
        normals.append(n)

    p = lat_lon_rad(-90,0,radius)
    n = scale3(p,1.0/mag(p))
    points.append(p)
    normals.append(n)
    return list(map( lambda x: add(x,center),points)),normals

# face indices for icosahedron
icaIndices = [ [1,11,3],[3,11,5],[5,11,7],[7,11,9],[9,11,1],
               [2,1,3],[2,3,4],[4,3,5],[4,5,6],[6,5,7],[6,7,8],[8,7,9],[8,9,10],[10,9,1],[10,1,2],
               [0,2,4],[0,4,6],[0,6,8],[0,8,10],[0,10,2] ]


def addVertex(nv,nn,verts,normals):
    for i in range(len(verts)):
        if vclose(nv,verts[i]):
            return i,verts,normals
    verts.append(nv)
    normals.append(nn)
    return len(verts)-1,verts,normals

def subdivide(f,verts,normals,rad):
    ind1 = f[0]
    ind2 = f[1]
    ind3 = f[2]
    v1 = verts[ind1]
    v2 = verts[ind2]
    v3 = verts[ind3]
    n1 = normals[ind1]
    n2 = normals[ind2]
    n3 = normals[ind3]
    va = add(v1,v2)
    vb = add(v2,v3)
    vc = add(v3,v1)
    ma = rad/mag(va)
    mb = rad/mag(vb)
    mc = rad/mag(vc)
    va = scale3(va,ma)
    vb = scale3(vb,mb)
    vc = scale3(vc,mc)
    
    na = add(n1,n2)
    na = scale3(na,1.0/mag(na))
    nb = add(n2,n3)
    nb = scale3(nb,1.0/mag(nb))
    nc = add(n3,n1)
    nc = scale3(nc,1.0/mag(nc))

    inda,verts,normals = addVertex(va,na,verts,normals)
    indb,verts,normals = addVertex(vb,nb,verts,normals)
    indc,verts,normals = addVertex(vc,nc,verts,normals)

    f1 = [ind1,inda,indc]
    f2 = [inda,ind2,indb]
    f3 = [indb,ind3,indc]
    f4 = [inda,indb,indc]

    return verts,normals, [f1,f2,f3,f4]

    

# make the sphere, return a surface representation
def sphere(diameter,center=point(0,0,0),depth=2):
    rad = diameter/2
    verts,normals = makeIcoPoints(center,rad)
    faces = icaIndices
    
    for i in range(depth):
        ff = []
        for f in faces:
            verts, norms, newfaces = subdivide(f,verts,normals,rad)
            ff+=newfaces
        faces = ff
                      
    
    return [verts,normals,faces]
    

        
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
    dd.linecolor = 'white'

    verts,normals,faces = sphere(50.0,point(0,0,0),3)
    
    for f in faces:
        dd.draw(line(verts[f[0]],
                     verts[f[1]]))
        dd.draw(line(verts[f[1]],
                     verts[f[2]]))
        dd.draw(line(verts[f[2]],
                     verts[f[0]]))

    dd.draw_surface(verts,normals,faces)
    dd.linecolor = 'aqua'
    dd.polystyle = 'points'
    dd.pointstyle = 'xo'
    dd.draw(verts)
    dd.display()
