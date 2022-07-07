# Changelog
This document tracks the notable changes to the LEOD project.

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