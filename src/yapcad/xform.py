## generalized matrix transformation operations for 3D homogeneous
## coordinates in yapCAD

## Copyright (c) 2020 Richard W. DeVaul
## Copyright (c) 2020 yapCAD contributors
## All rights reserved

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from math import *
import yapcad.geom as geom

## a matrix is represented as a list of four four vectors. In a
## matrix, vectors represent rows unless the transpose property is
## true.  Because vectors are represented as lists (not as instances
## of a class with meta-info) we assume that operations like Mx imply
## a column vector and that xM imply a row vector.

## There are two classes defined here: Matrix and MatrixStack.  Matrix
## is a relatively lightweight class that provides foundation
## matrix-matrix and matrix-vector operations.  MatrixStack is a class
## that captures a series of transformations and their inverses.  Note
## that inverses are generally determined analytically through
## composition, though we don't dismiss the possibility of a numeric
## inversion.

## FIXME: only one class is defined for now


class Matrix:
    """4x4 transformation matrix class for transforming homogemenous 3D coordinates"""

    def __init__(self,a=False,trans=False):
        self.m = [[1,0,0,0],
                  [0,1,0,0],
                  [0,0,1,0],
                  [0,0,0,1]]
        self.trans=False
        
        if isinstance(a,Matrix):
            for i in range(4):
                self.setrow(i,a.getrow(i))
                
        elif isinstance(a,(tuple,list)):
            if len(a) == 4:
                r1 =a[0]
                r2 =a[1]
                r3 =a[2]
                r4 =a[3]
                if len(r1) == len(r2) == len(r3) == len(r4) == 4:
                    for i in range(4):
                        for j in range(4):
                            x =a[i][j]
                            if (not isinstance(x,bool)) and \
                               isinstance(x,(int,float)):
                                self.m[i][j]=x
                            else:
                                raise ValueError('bad element in matrix initialization: {}'.format(x))
            elif len(a)==16:
                for i in range(4):
                    for j in range(4):
                        ind=i*4+j
                        x = a[ind]
                        if (not isinstance(x,bool)) and \
                           isinstance(x,(int,float)):
                            self.m[i][j]=x
                        else:
                            raise ValueError('bad element in matrix initialization: {}'.format(x))

            elif a == False:
                pass
            else:
                raise ValueError('bad thing used in attempt to initialize matrix: ()'.format(a))
        self.trans=trans
                            
    def __repr__(self):
        return "Matrix({},{},{},{},{})".format(self.m[0],self.m[1],
                                               self.m[2],self.m[3],self.trans)

    #return value indexed by i,j
    def get(self,i,j):
        if i < 0 or i > 3 or j < 0 or j > 3:
            raise ValueError('bad index passed to get: {},{}'.format(i,j)) 
        if self.trans:
            return self.m[j][i]
        else:
            return self.m[i][j]

    #set value indexed by i,j
    def set(self,i,j,x):
        if i < 0 or i > 3 or j < 0 or j > 3:
            raise ValueError('bad index passed to set: {},{}'.format(i,j))
        if geom.isgoodnum(x):
            if self.trans:
                self.m[j][i]=x
            else:
                self.m[i][j]=x
        else:
            raise ValueError('bad value passed to set: {}'.format(x))

    def getrow(self,i):
        if i < 0 or i > 3:
            raise ValueError('bad row passed to getrow: {}'.format(i))
        if self.trans:
            return [self.m[0][i],
                    self.m[1][i],
                    self.m[2][i],
                    self.m[3][i]]
        else:
            return self.m[i]
        
    def getcol(self,j):
        if j < 0 or j > 3:
            raise ValueError('bad column passed to getcol: {}'.format(j))
        if not self.trans:
            return [self.m[0][j],
                    self.m[1][j],
                    self.m[2][j],
                    self.m[3][j]]
        else:
            return self.m[j]

    def setrow(self,i,x):
        if not geom.isvect(x):
            raise ValueError('bad non-vector passed to setrow: {}'.format(x))
        if i < 0 or i > 3:
            raise ValueError('bad row index passed to setrow: {}'.format(i))
        if self.trans:
            self.m[0][i] = x[0]
            self.m[1][i] = x[1]
            self.m[2][i] = x[2]
            self.m[3][i] = x[3]
        else:
            self.m[i] = x
        
    def setcol(self,j,x):
        if not geom.isvect(x):
            raise ValueError('bad non-vector passed to setcol: {}'.format(x))
        if j < 0 or j > 3:
            raise ValueError('bad column index passed to setcol: {}'.format(j))
        if not self.trans:
            self.m[0][j] = x[0]
            self.m[1][j] = x[1]
            self.m[2][j] = x[2]
            self.m[3][j] = x[3]
        else:
            self.m[j] = x
        

    # matrix multiply.  If x is a matrix, compute MX.  If X is a
    # vector, compute Mx. If x is a scalar, compute xM. If x isn't any
    # of these, return False.  Respects transpose flag.
    
    def mul(self,x):
        if isinstance(x,Matrix):
            result = Matrix()
            for i in range(4):
                for j in range(4):
                    result.set(i,j,
                               geom.dot4(self.getrow(i),x.getcol(j)))
            return result
        elif geom.isvect(x):
            result = geom.vect()
            for i in range(4):
                result[i]=geom.dot4(self.getrow(i),x)
            return result
        elif geom.isgoodnum(x):
            result = Matrix()
            for i in range(4):
                result.setrow(i,geom.scale4(self.getrow(i),x))
            return result
        
        raise ValueError('bad thing passed to mul(): {}'.format(x))
    

# return the generalized 4x4 arbitrary axis rotation matrix
def Rotation(axis,angle,inverse=False):
    m = geom.mag(axis)
    u = axis
    if m < geom.epsilon:
        raise ValueError('zero-length rotation axis not allowed')
    if not geom.close(m,1.0):
        u = geom.geom.scale3(axis,1.0/m)

    if inverse:
        angle *= -1.0
    rad = (angle%360.0)*geom.pi2/360.0

    ux = u[0]
    uy = u[1]
    uz = u[2]
    
    cang = cos(rad)
    cmin = 1.0-cang
    
    sang = sin(rad)
    smin = 1.0-sang

    # # see https://en.wikipedia.org/wiki/Rotation_matrix
    ## THIS TURNS OUT TO BE WRONG! Bad wikipedia! no biscuit
    # R = [[cang + ux*ux*cmin, ux*uy*cmin-uz*sang, ux*uz*cmin+uy*sang,0],
    #      [uy*ux*cmin+uz*smin, cang + uy*uy*cmin, uy*uz*cmin - uz*smin,0],
    #      [uz*ux*cmin-uy*smin, uz*uy*cmin+ux*smin, cang+uz*uz*cmin,0],
    #      [0,0,0,1]]

    # see http://www.opengl-tutorial.org/assets/faq_quaternions/index.html#Q38
    ## This is correct
    R = [[cang + ux*ux*cmin, ux*uy*cmin-uz*sang, ux*uz*cmin+uy*sang,0],
         [uy*ux*cmin+uz*sang, cang + uy*uy*cmin, uy*uz*cmin - ux*sang,0],
         [uz*ux*cmin-uy*sang, uz*uy*cmin+ux*sang, cang+uz*uz*cmin,0],
         [0,0,0,1]]

    return Matrix(R)

def Translation(delta,inverse=False):
    if inverse:
        delta = mul(delta,-1.0)
    dx = delta[0]
    dy = delta[1]
    dz = delta[2]
    T = [[1,0,0,dx],
         [0,1,0,dy],
         [0,0,1,dz],
         [0,0,0,1]]
    return Matrix(T)

def Scale(x,y=False,z=False,inverse=False):
    sx = sy = sz = 1.0
    if geom.isgoodnum(x):
        sx = x
        if geom.isgoodnum(y) and geom.isgoodnum(z):
            sy = y
            sz = z
        else:
            sy = sz = x
    elif geom.isvect(x):
        sx = x[0]
        sy = x[1]
        sz = x[2]
    else:
        raise ValueError('bad scaling values passed to Scale')

    if inverse:
        sx = 1.0/sx
        sy = 1.0/sy
        sz = 1.0/sz
        
    S = [[sx,0,0,0],
         [0,sy,0,0],
         [0,0,sz,0],
         [0,0,0,1.0]]
    return Matrix(S)
    
