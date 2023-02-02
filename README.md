# ESD : Ellipsoid Shortest Distance {#mainpage}

## Overview {#overview}

ESD is a Python library for compuationally determining shortest distances over the surface of an ellipsoid.
This is primarily done using the fast marching method (FMM) of Kimmel and Sethian \cite Kimmel1998, on a suitably defined mesh 
of the ellipsoid surface. This approach allows the computation of approximate distances between any pair of points on 
any ellipsoid (sphere, spheroid or triaxial) with great efficiency. Moreover, multiple distances can be computed from 
a particlular start point using a single calculation. This makes the method particularly well suited to solving single source problems.

To use this library, one must first ensure that the ESD/ directory is added to the PYTHONPATH. One can then import and use the 
code implemented in the esd/ directory within Python programs.

This library uses both NumPy and SciPy (as well as Matplotlib in the data/figures/ directory). The GeographicLib library of 
Karney \cite Karney2011 \cite Karney2012 \cite Karney2019 is also used for making comparisons, but is not required for the majority of the program.

The examples/ directory contains a set of Python scripts that illustrate how the ESD library is used.

To read the full documentation, open index.html in the docs/html/ directory.

## Features {#features}

The ESD library includes code implementing the following:
- Generation of surface meshes on the ellipsoid surface:
    - Using polar coordinates to define a parametric mesh.
	- Using a regular icosahedron to define a triangulation of the ellipsoid surface.
- Dijkstra's algorithm and the fast marching method on an ellipsoid surface mesh:
    - Between a pair of points.
	- From a given start point to multiple end points.
	- From a given start point to all mesh vertices, allowing interpolation of distance to any end point after the initial calculation.
- Conversions between Cartesian and polar coordinates on the surface of a given ellipsoid.
- Conversions between Cartesian and ellipsoidal coordinates on the surface of a given ellipsoid.
- Finding the ellipse of intersection between a plane and an ellipsoid.
- Alternative algorithms for computing distances on the ellipsoid:
    - True distance via an analytical formula (only possible for spheres).
	- Spheroid geodesic distances using GeographicLib.
	- Spheroid and triaxial geodesic paths/distances using the boundary value method (BVM) of Panou et al. \cite Panou2013a \cite Panou2013.
	- Taxicab distances (i.e. rectilinear distance with respect to polar coordinates).

## Author {#author}

   * Callum Marples