### User-configurable parameters (some of these should be moved to being command line inputs)
spw = '110_spw29' 
prefix = '/blue/adamginsburg/abulatek/brick/first_results/temperature_map/' # For temperature and column density maps

# Import various packages
import numpy as np
from spectral_cube import SpectralCube 
from astropy import units as u
from astropy.io import fits
import pyspeckit
import pylab as pl
from pyspeckit.spectrum.models import lte_molecule

# Retrieve the data
import requests
# 10 MB file: stream it, but it should be able to load in memory
response = requests.get('https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/8QJT3K/POCX2W', stream=True)
with open('temp_cube.fits', 'wb') as fh:
    fh.write(response.content)

# Set up things based on user-input parameters
cube_filename = f'/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/source_ab_{spw}_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(cube_filename, format='casa_image')
print(cube)
freq_axis = cube.spectral_axis

# Import temperature and column density maps and generate moment 1 map (for inputting fitted values)
temp = fits.getdata(f'{prefix}temperature_map_test.fits')
log_N = fits.getdata(f'{prefix}log_N_test.fits')
masked_cube = SpectralCube.read(f'{prefix}methyl_cyanide/ch3cn_0_masked.fits', format='fits')
mom1 = masked_cube.moment1()

# Take a spectrum from the cube
meanspec = cube.max(axis=(1,2)) #, how='slice', progressbar=True)
meanspec_K = meanspec.value * cube.jtok_factors() # have to hack this, which is _awful_ and we need to fix it
sp = pyspeckit.Spectrum(xarr=cube.spectral_axis, data=meanspec_K)

# The LTE model doesn't include an explicit filling factor, so we impose one
fillingfactor = 0.05

# There's continuum in our spectrum, but it's also not modeled
offset = np.nanpercentile(sp.data.filled(np.nan), 20)

# Look at a couple species
# Be cautious! Not all line catalogs have all of these species, some species
# can result in multiple matches from splatalogue, and definitely not all
# species have the same column and excitation
species_list = ('HNCO', 'CH3CN',) # 'C18O', '13CO', #, 'HNCO'
# REMEMBER THE COMMA AT THE END TO MAKE THIS A TUPLE!
# Molecules that didn't work: 'H13CN', 'H13CO+', 'HN13C', 'C18O', '13CO',
# I think I don't know how to search for these, and I don't know which database I should be picking names from. Need to find that out

# Make a set of subplots
# Do we need something else here? The first time this gets run in a notebook, nothing plots
fig, axes = pl.subplots(len(species_list)+1, 1, figsize=(18,14), num=1)
pl.draw()

# Generate model spectra for each molecule individually, then add them all up (also plot all of this)
mods = []
for species, axis in zip(species_list, axes):
    print(species, axis)
    freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters_JPL(species, 
                                                                              fmin=sp.xarr.min(),
                                                                              fmax=sp.xarr.max())
    fine_xarr = np.linspace(sp.xarr.min(), sp.xarr.max(), 1000)
    mod = lte_molecule.generate_model(sp.xarr, # In order to incr. spectral res, need this to be fine
                                      50*u.km/u.s, 
                                      1.5*u.km/u.s, 
                                      100*u.K, 
                                      5e15*u.cm**-2, 
                                      freqs, aij, deg, EU, partfunc) # Can add another velocity component
    mods.append(mod)

    sp.plotter(axis=axis, figure=fig, clear=False)
    sp.plotter.axis.plot(sp.xarr, mod*fillingfactor+offset, label=species, zorder=-1)
    sp.plotter.axis.text(0.5,0.9,species,transform=axis.transAxes)

    if axis != axes[-1]:
        axis.set_xlabel("")
        axis.set_xticklabels([])
    if axis != axes[len(species_list) // 2]:
        axis.set_ylabel("")

axis = axes[-1]
sp.plotter(axis=axis, figure=fig, clear=False)
sp.plotter.axis.plot(sp.xarr, np.nansum(mods,axis=0)*fillingfactor+offset, label='Sum', zorder=-1) #orig sp.xarr
sp.plotter.axis.text(0.5,0.9,'Sum',transform=axis.transAxes)
