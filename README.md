# Surface-surface intersection curves

The code is based on the method of Bo, Liu and Zhang, Topological
structure extraction for computing surface–surface intersection
curves, The Visual Computer, 2025,
https://doi.org/10.1007/s00371-024-03616-8 and the crawl algorithm of
Parakkat and Muthuganapathy, Crawl through neighbors: a simple curve
reconstruction algorithm, Eurographics Symposium on Geometry
Processing 2016, https://dx.doi.org/10.1111/cgf.12974

# Code from elsewhere

This code includes modified versions of
* Ramsey, Potter, and Hansen's ray-bilinear path intersection code
  from http://shaunramsey.com/research/bp/, originally released under
  the MIT License. Modification is translation to C.
* John Burkardt's C version of Richard Brent's PRAXIS code for
  Powell's method, from
  https://people.math.sc.edu/Burkardt/c_src/praxis/praxis.html,
  released under the LGPL. The modification is the addition of
  parameter passing to the objective function, and some reformatting.

The Delaunay triangulation code uses Jonathan Shewchuk's Triangle code
in the library version of Christian Woltering,
https://github.com/wo80/Triangle, which you need to install before
compiling the code. 

# Installation

There is a simple Makefile which should be sufficient to compile the
code

`make`

If you have installed libraries or include files somewhere not
normally searched you may need to specify their location when you call
make. For example

`CFLAGS=-I[location of include files] LDFLAGS=-L[location of libraries] make`

# Running the test code

To run the test code

`./ssi-test`

This will generate two surfaces and calculate their intersections. The
surfaces are output to files `bvh1.dat` and `bvh2.dat`, points on the
intersection curve are output to `points.dat` and the point lists are
output to `curves.dat`.

Surfaces are output as a list of points on a square grid with each
line giving

`x y z u v`

where the $(x,y,z)$ is the spatial position of a node and the surface
parametric coordinates are $(u,v)$.

Intersection points are output as 

`x1 y1 z1 x2 y2 z2 u1 v1 u2 v2`

where $(x_i,y_i,z_i)$ is position calculated in surface $i$ and the
parametric coordinates $(u_i,v_i)$ are those in the $(0,1)\times(0,1)$
space. 
