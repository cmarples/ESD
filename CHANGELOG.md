# Changelog
This document tracks the notable changes to the LEOD project.

## 2022-11-25
### Added
- Started on script to plot the resolution date (figure_resolution.py)
### Modified
- Bug in Dijkstra's algorithm fixed (not adding distance to trial point).

## [0.0.3] - 2022-11-24
### Added
- Accuracy as a function of resolution started.
- Alternate methods applied for resolution comparison.
### Modified
- Attempts at refinement and 2nd order FMM abandoned due to complexity and time required for only marginal accuracy gains.
- Precalculation routine no longer outputs a maximum angle.
- Polar grid generation routine no longer takes a 'Dijkstra method flag' as an input (as unnecessary).

## 2022-11-17
### Added
- Edge splitting on the icosahedral triangulation implemented.
- One-sided neighbours caused by these splits cause a bug in the FMM code that needs fixing.

## 2022-11-16
### Added
- Plot of icosahedral triangulation, with face colours according to the largest angle.

## 2022-11-15
### Modified
- Icosahedron based triangulation updated for compatibility with newest FMM code.
- Bug in bottom left vertex of triangular faces fixed.

## 2022-11-14
### Added
- Beginnings of a script to compute geodesics on a set of test examples.
- Finding obtuse angles relevant to a given pair distance.
### Modified
- Single endpoint fast marching calculation now performed within separate routine.

## 2022-11-11
### Added
- Early termination of FMM if all endpoints reached (closest vertex and neighbours accepted).

## 2022-11-10
### Added
- Interpolation of endpoint distance, using inverse distance weights.

## 2022-11-9
### Modified
- Optimisation of FMM by precalculating all grid data (neighbour distances, faces and face angles)

## [0.0.3] - 2022-11-8
### Modified
- Triangulation approach likely to be long, both in implementation and in run-time.
  Thus the theta-phi grid has been returned in the hope that obtuse angles have little effect.
- Started optimisation of the fats marching code.


## 2022-10-31
### Added
- Routine for finding closest vertex on a triangulation to a given point (using a descending direction approach)

## 2022-10-29
### Modified
- Sphere triangulation completed.

## 2022-10-28
### Added
- Started implementation of FMM on a general triangulation.
### Modified
- Removed some old code (placed in seperate directory).

## 2022-10-25
### Added
- Attempt at obtuse angle splitting on triaxial ellipsoid (theta-phi grid).

## 2022-10-18
### Added
- First order fast marching update in two forms: 'trigonometric form' and 'linear algebra' form.
- Both forms debugged and giving the same answers.

## 2022-10-10
### Added
- Function for finding distance between neighbours.

## 2022-10-09
### Modified
- Improved face finding procedure.

## 2022-10-08
### Added
- Face finding for 4 or 8 neighbours.

## 2022-10-05
### Modified
- Grid class removed. Using a list of 'vertex' objects instead.

## [0.0.2] - 2022-10-04
### Modified
- FMM code restarted to allow use on general triangulations as well as on a theta-phi grid.

## 2022-10-03
### Added
- Triangulation of a sphere (using faces of an icosahedron) implemented.

## 2022-09-21
### Added
- Multiple endpoint FMM implemented.

## 2022-09-13
### Added
- Dijkstra's Algorithm with 8 neighbours implemented.

## 2022-09-12
### Added
- Boundary Value Method (of Panou 2013) implemented in LEOD for triaxial geodesics.

## 2022-09-5
### Added
- GeographicLib calls for spheroid geodesics implemented within LEOD.

## 2022-08-30
### Added
- Test example run with varying grid aspect ratio.

## 2022-08-25
### Added
- Simple FMM on a flat grid implemented for testing purposes.

## 2022-08-22
### Modified
- Conditions for grid downsampling revised.

## 2022-08-19
### Modified
- Additional bugs in source refinement procedure fixed.

## 2022-08-18
### Modified
- Bug in border pixel determination fixed.
- Bug with refinement at low resolution fixed.

## 2022-08-14
### Added
- Flags used to determine if a pixel is on the refined region border.
### Modified
- 'pixel' reimplemented as a dictioary (so only those needed for refinement are defined).

## 2022-08-13
### Modified
- 'alive' and 'geo_distances' reimplemented as dictionaries for speed and simplicity.

## 2022-08-07
### Modified
- Bug in theta finding fixed.
- Incorrect handling of 'alive' in transfer fixed.

## 2022-08-02
### Modified
- Refined pixel finding sped up.

## 2022-07-28
### Modified
- Neighbour search (when creating refined region) sped up.

## 2022-07-27
### Added
- Source refined FMM implementation.

## 2022-07-22
### Added
- Mapping from refined grid to main grid (transfer_grid).
- Calculation of subsequent main grid fast marching method.

## 2022-07-20
### Added
- Neighbours of refined pixels (in the refinement region) now initialised.
### Modified
- Fixed potential bug where south pole pixel incorrectly indexed as n_theta-1, rather than n_pixels-1 when adding to 'border_pixels'.
- Renamed 'neighbour_distance' function to 'get_neighbour_distance' to avoid ambiguity with the GeoPixel lists of the same name.

## 2022-07-16
### Added
- Border list in refined grid has been implemented.

## 2022-07-15
### Modified
- GeoFMM functions now take a GeoGrid as an input, so that a refined grid can be used.

## 2022-07-14
### Added
- First part of initialisation of a source refined grid.

## 2022-07-11
### Added
- Implementation of the 4-neighbour Dijkstra algorithm (started 2022-07-08).
- Implementation of the first order fast marching method.
- Implementation of the second order fast marching method.
- Testing of 4-neighbour Dijkstra, 1st and 2nd order FMM; comparing to test from EORL.

## 2022-07-07
### Added
- Determination of neighbouring pixels
- Function to get (or calculate) distance between neighbouring pixel centres.

## 2022-07-06
### Modified
- Changed the name of the project to 'LEOD', Library for Ellipsoid cOmputational Determination.
- Modified the directory structure, with the contents of the 'leod' subdirectory as a package.
### Added
- apps directory containing scripts to be executed (main.py renamed example_geodesics.py).
- tests directory containing unit tests (test.py renamed test_geodesics.py).

## 2022-07-05
### Added
- geo_grid.py : Grid initialisation and index finder routines.
- geo_pixel.py : Pixel initialisation.
- main.py : Script intended to test run the geodesic routines.
- test.py : Unit testing script (including index finder tests).

## [0.0.1] - 2022-07-04
### Added
- ellipsoid_shape.py, defining the ellipsoid_shape class.
- An empty __init__.py file.
- The README file.
- The AUTHORS file.
- This CHANGELOG file.