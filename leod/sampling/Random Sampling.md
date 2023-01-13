# Random Sampling of Surfaces

### Overview {#overview}

The sampling directory contains classes used to generate random points on the surface of a sphere (sphere.py) or an ellipsoid (ellipsoid.py).
Each class corresponds to a particular sampling algorithm. All of these classes have a common method called random_surface_point(), which returns
a single randomly generated point on the surface of the relevant shape.

### Sphere Samplers

The file sampling/sphere.py contains classes for the following six sphere surface sampling algorithms: 

* Cubic rejection.
* Marsaglia's method.
* Cook's method.
* Gaussian random points.
* The trigonometric method.
* Area rejection.

### Ellipsoid Samplers

The file sampling/ellipsoid.py contains classes for the following five ellipsoid surface sampling algorithms: 

* Naive scaling (does not provide uniformly random points).
* Gradient rejection.
* Area rejection (polar coordinates).
* Area rejection (Mercator coordinates).
* The 'generic' sampler of Detwiler et al.

### Sources of Random Numbers

To generate uniform random numbers, two algorithms from the NumPy library were included in the sampler classes:

* Mersenne twister (MT19937).
* Permuted Congruential Generator (PCG64).

By allowing a choice between more than one source of random numbers, it is possible to run tests with each to check that the distribution
of sampled sphere/ellipsoid surface points is not influenced by the random number generator.










