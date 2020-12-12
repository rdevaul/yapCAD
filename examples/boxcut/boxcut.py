## boxcut procedural design example for yapCAD
## Born on 12 December, 2020
## Copyright (c) 2020 Richard DeVaul

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

"""
procedural design example that generates squeeze-fit box designs from
input parameters and produces DXF files as output.
"""

from yapcad.geom import *
from yapcad.poly import *
from yapcad.combine import *

if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description="""
    procedural design example that generates squeeze-fit box designs from
    input parameters and produces DXF files as output.
    """)

    parser.add_argument("length",type=float, help="length of the box")
    parser.add_argument("width",type=float, help="width of the box")
    parser.add_argument("height",type=float, help="height of the box")
    parser.add_argument("-u","--units",type=str,choices=["mm","inch"],
                        default="mm",help="units for dimensions")
    parser.add_argument("-o","--outname",type=str,default="boxout",
                        help="file name base for output, default is boxout")
    parser.add_argument("-s","--separate",action="store_true",
                        help="make one DXF drawing per face rather than a combined tiled layout")
    parser.add_argument("-t","--thick",type=float,default=-1,
                        help="material thickness, default is 3.175 mm (1/4\")")
    parser.add_argument("-k","--kerf",type=float,default=-1,
                        help="kerf compensation, default is 0.397mm (1/64\")")
    parser.add_argument("-g","--gl",action="store_true",
                        help="do an interactive 3D visualization of the box")
    args = parser.parse_args()
    print(args)
    print("done")

