from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astroquery.xmatch import XMatch

Vizier.ROW_LIMIT = -1

# Retrieve VVV catalog in a radius around the Brick
VVV_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'),
                                radius=1*u.arcmin,
                                catalog='II/348/vvv2')[0]

print("The VVV catalog around the Brick has this many sources:", len(VVV_B))

# Cross-match the GALACTICNUCLEUS survey with Gaia DR2
VVV_GDR2_xmatch = XMatch.query(cat1=VVV_B, cat2='vizier:I/345/gaia2', max_distance=1*u.arcsec, colRA1='RAJ2000', colDec1='DEJ2000', colRA2='RA_ICRS', colDec2='DE_ICRS')
print("The VVV/GDR2 crossmatch catalog has this many sources:", len(VVV_GDR2_xmatch))
print("That means there are this many sources in VVV that aren't in GDR2:", len(VVV_B) - len(GN_GDR2_xmatch))
