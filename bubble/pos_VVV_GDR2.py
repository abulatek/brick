from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import sys

# Import the catalog from the VO table file
catalog = Table.read('OB_star_candidates_shift.xml')

# Access the coordinates for the stars, already in decimal degrees
RA = catalog['RAJ2000']
DE = catalog['DEJ2000']

# Convert from RA/DE to galactic coordinates
coords = SkyCoord(ra=RA*u.degree, dec=DE*u.degree, frame='icrs')
coords_galactic = coords.galactic
l = coords_galactic.l
b = coords_galactic.b

# Plot position of sources on the sky
plt.figure(dpi = 100)
plt.grid()
plt.scatter(RA, DE, color='xkcd:blue', alpha=0.5)
plt.xlabel("RA (degrees)")
plt.ylabel("DE (degrees)")
plt.title("Locations of candidate OB stars from VVV DR2", size = 14)
plt.show()

# Plot position of sources on the sky
plt.figure(dpi = 100)
plt.grid()
plt.scatter(l, b, color='xkcd:blue', alpha=0.5)
plt.xlabel("l (degrees)")
plt.ylabel("b (degrees)")
plt.title("Locations of candidate OB stars from VVV DR2", size = 14)
plt.show()
