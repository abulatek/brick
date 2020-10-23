from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astroquery.xmatch import XMatch

Vizier.ROW_LIMIT = -1

# Retrieve GALACTICNUCLEUS central catalog in the 1 arcminute radius region around the Brick
GN_central_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'), 
                                radius=10*u.arcmin, catalog='J/A+A/631/A20/central')[0]

# Retrieve GALACTICNUCLEUS nsdeast catalog in the 1 arcminute radius region around the Brick
GN_nsdeast_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'), 
                                radius=10*u.arcmin, catalog='J/A+A/631/A20/nsdeast')[0]

'''
# Try the other five sub-catalogs
GN_nsdwest_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=10*u.arcmin,
                                catalog='J/A+A/631/A20/nsdwest')[0]
GN_teast_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=10*u.arcmin,
                                catalog='J/A+A/631/A20/teast')[0]
GN_twest_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=10*u.arcmin,
                                catalog='J/A+A/631/A20/twest')[0]
GN_ibnorth_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=10*u.arcmin,
                                catalog='J/A+A/631/A20/ibnorth')[0]
GN_ibsouth_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=10*u.arcmin,
                                catalog='J/A+A/631/A20/ibsouth')[0]
'''

# Stack the sources in the GN catalogs together
GN_B_stack = vstack([GN_central_B, GN_nsdeast_B])

print("The GN central catalog 1 arcmin around Brick has this many sources:", len(GN_central_B))
print("The GN nsdeast catalog 1 arcmin around Brick has this many sources:", len(GN_nsdeast_B))
print("The stacked GN catalog around the Brick has this many sources:", len(GN_B_stack))

'''
#print("nsdwest:", len(GN_nsdwest_B))
#print("teast:", len(GN_teast_B))
#print("twest:", len(GN_twest_B))
#print("ibnorth:", len(GN_ibnorth_B))
#print("ibsouth:", len(GN_ibsouth_B))
'''

# Cross-match the GALACTICNUCLEUS survey with Gaia DR2

# may need to convert ra/dec into decimal degrees

GN_GDR2_xmatch = XMatch.query(cat1=GN_B_stack, cat2='vizier:I/345/gaia2', max_distance=1*u.arcsec, colRA1='RAJ2000', colDec1='DEJ2000', colRA2='RA_ICRS', colDec2='DE_ICRS')
print("The GN/GDR2 crossmatch catalog has this many sources:", len(GN_GDR2_xmatch))
print("That means there are this many sources in GN that aren't in GDR2:", len(GN_B_stack) - len(GN_GDR2_xmatch))
