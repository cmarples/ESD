print(f'Invoking __init__.py for {__name__}')
import leod.ellipsoid_shape
import leod.fmm_polar_graph
import leod.triangulation_sphere
import leod.fmm_precalculation
import leod.fmm_callers
import leod.sphere_geodesics
import leod.triaxial_geodesics
import leod.taxicab_distance
# Import GeographicLib, if it has been installed
try:
    from geographiclib.geodesic import Geodesic
    import leod.spheroid_geodesics
except:
    print("GeographicLib not installed.")