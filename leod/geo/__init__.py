import leod.geo.sphere
import leod.geo.triaxial
import leod.geo.taxicab
# Import GeographicLib, if it has been installed
try:
    from geographiclib.geodesic import Geodesic
    import leod.geo.spheroid
except:
    print("GeographicLib not found.")