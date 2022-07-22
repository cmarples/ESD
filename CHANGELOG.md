# Changelog
This document tracks the notable changes to the LEOD project.

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