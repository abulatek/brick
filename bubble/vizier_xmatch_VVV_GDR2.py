from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, setdiff
from astroquery.xmatch import XMatch
import sys

Vizier.ROW_LIMIT = -1

# Retrieve VVV catalog in a radius around the Brick
VVV_B = Vizier.query_region(SkyCoord(ra=266.544, dec=-28.7049, unit=(u.deg, u.deg), frame='icrs'), 
                                radius=1*u.arcmin, 
                                catalog='II/348/vvv2')[0]

print("The VVV catalog around the Brick has this many sources:", len(VVV_B))

# Cross-match the GALACTICNUCLEUS survey with Gaia DR2
VVV_GDR2_xmatch = XMatch.query(cat1=VVV_B, cat2='vizier:I/345/gaia2', max_distance=1*u.arcsec, colRA1='RAJ2000', colDec1='DEJ2000', colRA2='RA_ICRS', colDec2='DE_ICRS')
print("The VVV/GDR2 crossmatch catalog has this many sources:", len(VVV_GDR2_xmatch))
print("That means there are this many sources in VVV that aren't in GDR2:", len(VVV_B) - len(VVV_GDR2_xmatch))

# Keep only the stars in VVV that are NOT in the crossmatch, using the column 'iauname' as the link between the catalogs
sdiff = setdiff(VVV_B, VVV_GDR2_xmatch, keys=['iauname'])
#print(len(sdiff))
#print(sdiff)
#print(sdiff.colnames)

# Suggestion from Adam on how to remove masked values
#print(len(sdiff['Ksmag3'].mask))

# Send the working catalog to pandas to be able to drop rows with NaN values in one or more columns
sdiff_df = sdiff.to_pandas()

# Keep only stars that have measured K-band magnitudes
sdiff_df_dropped = sdiff_df.dropna(subset=['Ksmag3'])
print("This many sources are detected in the Ks band:", len(sdiff_df_dropped))

# Move back to astropy tables
sdiff_dropped = Table.from_pandas(sdiff_df_dropped)

# Filter for sources that are too bright
# K band cutoff:
K_mask = sdiff_dropped['Ksmag3'] > 7.0
sdiff_dropped_K = sdiff_dropped[K_mask]
print("This many sources are appropriately bright in the Ks band:", len(sdiff_dropped_K))

# J band cutoff:
J_mask = sdiff_dropped_K['Jmag3'] > 9.76
sdiff_dropped_K_J = sdiff_dropped_K[J_mask]
print("This many sources are appropriately bright in the J band:", len(sdiff_dropped_K_J))

# H band cutoff:
H_mask = sdiff_dropped_K_J['Hmag3'] > 8.13
sdiff_dropped_K_J_H = sdiff_dropped_K_J[H_mask]
print("This many sources are appropriately bright in the H band:", len(sdiff_dropped_K_J_H))
