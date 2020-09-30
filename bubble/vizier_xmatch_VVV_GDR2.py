from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, setdiff
from astroquery.xmatch import XMatch
import sys

Vizier.ROW_LIMIT = -1

# Retrieve VVV catalog in a radius around the Brick
VVV_B = Vizier.query_region(SkyCoord(l=0.2509749, b=0.01616672, unit=(u.deg, u.deg), frame='galactic'), 
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

# Remove sources with *perrbits values that are nonzero for the J, H, and Ks filters
# This is not the best way to do this, but I will fix it in the future when I have more time
perrbits_mask_J = sdiff_dropped_K_J_H['Jperrbits'] == 0 
sdiff_reliable_J = sdiff_dropped_K_J_H[perrbits_mask_J]
perrbits_mask_H = sdiff_reliable_J['Hperrbits'] == 0
sdiff_reliable_J_H = sdiff_reliable_J[perrbits_mask_H]
perrbits_mask_Ks = sdiff_reliable_J_H['Ksperrbits'] == 0
sdiff_reliable_J_H_K = sdiff_reliable_J_H[perrbits_mask_Ks]
print("This many sources have *perrbits equal to 0:", len(sdiff_reliable_J_H_K))

# Export the intermediate catalog to a VO table file
#sdiff_reliable_J_H_K.write('OB_star_candidates_shift.xml', table_id='OB_star_catalog_shift', format='votable')

# Pick out the 16 brightest stars in the Ks band
sdiff_reliable_J_H_K.sort(['Ksmag3'])
catalog_subset = sdiff_reliable_J_H_K[:16] # First row is the header, I think
print("This many sources are in the final subset catalog:", len(catalog_subset))

# Export the subset catalog to a VO table file
#catalog_subset.write('OB_star_candidates_subset_shift.xml', table_id='OB_star_catalog_subset_shift', format='votable')