from astropy.table import Table
import matplotlib.pyplot as plt
import sys

# Import the catalog from the VO table file
catalog = Table.read('OB_star_candidates.xml')

# Access the photometric values individually; they already have units of mag
Jmag3 = catalog['Jmag3']
Hmag3 = catalog['Hmag3']
Ksmag3 = catalog['Ksmag3']

# Plot color-magnitude diagrams
plt.figure(dpi = 100)
plt.grid()
plt.scatter(Jmag3-Ksmag3, Jmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("J - Ks")
plt.ylabel("J")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

plt.figure(dpi = 100)
plt.grid()
plt.scatter(Jmag3-Ksmag3, Ksmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("J - Ks")
plt.ylabel("Ks")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

plt.figure(dpi = 100)
plt.grid()
plt.scatter(Jmag3-Hmag3, Jmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("J - H")
plt.ylabel("J")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

plt.figure(dpi = 100)
plt.grid()
plt.scatter(Jmag3-Hmag3, Hmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("J - H")
plt.ylabel("H")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

plt.figure(dpi = 100)
plt.grid()
plt.scatter(Hmag3-Ksmag3, Hmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("H - Ks")
plt.ylabel("H")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

plt.figure(dpi = 100)
plt.grid()
plt.scatter(Hmag3-Ksmag3, Ksmag3, color='xkcd:black', alpha=0.5)
plt.xlabel("H - Ks")
plt.ylabel("Ks")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.show()

