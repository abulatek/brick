from astropy.table import Table
import sys

catalog = Table.read('OB_star_candidates_subset.xml')
print("ORIGINAL COORDS")
print(catalog)

catalog_shift = Table.read('OB_star_candidates_subset_shift.xml')
print("SHIFTED COORDS")
print(catalog_shift)
