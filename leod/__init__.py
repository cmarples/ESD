print(f'Invoking __init__.py for {__name__}')
import leod.shape
import leod.fmm
import leod.geo


# Import GeographicLib, if it has been installed
try:
    from geographiclib.geodesic import Geodesic
    import leod.spheroid_geodesics
except:
    print("GeographicLib not installed.")