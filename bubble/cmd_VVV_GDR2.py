from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import sys

plt.rc('axes', axisbelow=True) # Put grid lines behind data points
plt.rc('text', usetex = True) # Use LaTeX font in plots
plt.rcParams['text.latex.preamble'] = [r'\usepackage{gensymb}']
#                                       r'\usepackage{sansmath}',
#                                       r'\sansmath']

# Import the catalog from the VO table file
catalog = Table.read('OB_star_candidates_shift.xml')

# Access the photometric values individually; they already have units of mag
Jmag3 = catalog['Jmag3']
Hmag3 = catalog['Hmag3']
Ksmag3 = catalog['Ksmag3']

# Plot color-magnitude diagrams 
'''
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
'''

# Grab the selected sources
catalog.sort(['Ksmag3'])
catalog_subset_0 = catalog[:16]
print(catalog_subset_0)
remove_mask_1 = catalog_subset_0['Ksmag3'] != 11.204
catalog_subset_1 = catalog_subset_0[remove_mask_1]
remove_mask_2 = catalog_subset_1['Ksmag3'] != 13.117
catalog_subset_2 = catalog_subset_1[remove_mask_2]
remove_mask_3 = catalog_subset_2['Ksmag3'] != 13.433
catalog_subset_3 = catalog_subset_2[remove_mask_3]

Jmag3_sel = catalog_subset_3['Jmag3']
Hmag3_sel = catalog_subset_3['Hmag3']
Ksmag3_sel = catalog_subset_3['Ksmag3']

plt.figure(dpi = 150)
plt.grid()
plt.scatter(Hmag3-Ksmag3, Ksmag3, color='xkcd:black', marker='+', label="catalog")
plt.scatter(Hmag3_sel-Ksmag3_sel, Ksmag3_sel, color='xkcd:azure', marker='o', label="to observe")
plt.gca().invert_yaxis()
plt.xlabel(r"$H - K_s$ (color)", size = 12)
plt.ylabel(r"$K_s$ (apparent magnitude)", size = 12)
plt.ylim(17, 11)
plt.legend(loc="upper left")
plt.title("Color-magnitude diagram for candidate OB stars from VVV DR2", size = 14)
plt.savefig("cmd.png")
plt.show()
