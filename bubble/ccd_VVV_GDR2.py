from astropy.table import Table
import matplotlib.pyplot as plt
import sys

# Import the catalog from the VO table file
catalog = Table.read('OB_star_candidates_shift.xml')

# Access the photometric values individually; they already have units of mag
Jmag3 = catalog['Jmag3']
Hmag3 = catalog['Hmag3']
Ksmag3 = catalog['Ksmag3']

# Plot color-color diagram
plt.figure(dpi = 100)
plt.grid()
plt.scatter(Hmag3-Ksmag3, Jmag3-Hmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("H - Ks")
plt.ylabel("J - H")
plt.title("Color-color diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()
