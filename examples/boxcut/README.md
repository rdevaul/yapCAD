# boxcut example
An example procedural design system in **yapCAD** for creating
squeeze-fit box designs cut from flat sheets of material

## how to use

install **yapCAD**, and invoke the design system as follows:

	python3 boxcut <length mm> <width mm> <height mm> [options]
	
where options includes:

	outname=<dxf file name base> # default is boxout
	allinone={true,false} # make one combined DXF file (default), or a 
	                      # separate file for each piece
	thick=<sheet thickness mm> # default is 3.175 mm, or 0.125"
	kerf=<kerf mm> # default is 0.397, or 1/64"
	opengl={true,false} # default is false -- visaulize result in 3D
	
	
For example, let's say we want to create a box that is 6" in length,
4" in width, and 3" in height. We are cutting quarter-inch stock
acrylic (the default for `boxcut`) on a laser cutter without built-in
kerf compensation, so we need to specify the kerf, which happens to be
about 0.5 mm.

The first thing to do is to work out the dimensions of the box in mm,
by multiplying inches by 25.4. That works out to 152.4mm x 101.6mm x
76.2mm. Next, let's decide what name to give the output design file
&mdash; let's call this design boxout6x4x3. Finally, we want to see
the result visualized in 3D, so we will need to select that option.

Ok, we are ready run boxcut as follows:

	python3 boxcut 252.4 101.6 76.2 outname=boxout6x4x3 kerf=0.5 opengl=true
	
This will cause an OpenGL window to pop up with a visualization of the
box, including instructions for manipulating the view.  Hit `esc` to exit.

The design file boxout6x4x3-combined.dxf will also be produced.  When
you open this file in your favorite DXF file viewer, you will see the
profile for each of the six pieces tiled and labeled.

The top surface of each piece as cut represent the outside surface of
each piece as assembled

