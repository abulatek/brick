import numpy as np

# retrieve the data
import requests

# 10 MB file: stream it, but it should be able to load in memory
response = requests.get('https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/8QJT3K/POCX2W', stream=True)

with open('temp_cube.fits', 'wb') as fh:
    fh.write(response.content)

from spectral_cube import SpectralCube 
from astropy import units as u
# there is some junk at low frequency
cube = SpectralCube.read('temp_cube.fits').spectral_slab(218.14*u.GHz, 218.54*u.GHz)

# this isn't really right b/c of beam issues, but it's Good Enough
meanspec = cube.mean(axis=(1,2))
# have to hack this, which is _awful_ and we need to fix it
meanspec_K = meanspec.value * cube.beam.jtok(cube.spectral_axis)

import pyspeckit

sp = pyspeckit.Spectrum(xarr=cube.spectral_axis, data=meanspec_K)


# alyssa is testing
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/spectra/'
# freq_spw = '139_spw71'
freq_spw = '146_spw51'
max_fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
# mean_fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
from spectral_cube import OneDSpectrum
from astropy.io import fits
kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn)).to(u.K)
sp = pyspeckit.Spectrum.from_hdu(kspectrum.hdu)
sp.xarr.convert_to_unit('GHz')


# The LTE model doesn't include an explicit filling factor, so we impose one
fillingfactor = 0.05

# there's continuum in our spectrum, but it's also not modeled
offset = np.nanpercentile(sp.data.filled(np.nan), 20)

# look at a couple species
# Be cautious! Not all line catalogs have all of these species, some species
# can result in multiple matches from splatalogue, and definitely not all
# species have the same column and excitation
# species_list = ('CH3OH', 'HCCCN', 'H2CO', 'CH3OCHO', 'CH3OCH3', 'C2H5CN', 'CN, v = 0, 1', 'H2CS',) #CNv=0?
# species_list = ('H2CS',)
species_list = ('CH3CN',)

# make a set of subplots:
import pylab as pl
# pl.close(1)
fig, axes = pl.subplots(len(species_list)+1, 1, figsize=(18,14), num=1)
pl.draw()

from pyspeckit.spectrum.models import lte_molecule
mods = []
for species, axis in zip(species_list, axes):
    freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters_JPL(species, fmin=sp.xarr.min(), fmax=sp.xarr.max())
    #freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(species, fmin=200*u.GHz, fmax=250*u.GHz, export_limit=1e5, molecule_name_jpl=species_jpl,
    #                                                                      line_lists=['SLAIM'])
    mod = lte_molecule.generate_model(sp.xarr, 50*u.km/u.s, 3*u.km/u.s, 100*u.K, 1e17*u.cm**-2, freqs, aij, deg, EU, partfunc)
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
sp.plotter.axis.plot(sp.xarr, np.nansum(mods,axis=0)*fillingfactor+offset, label='Sum', zorder=-1)
sp.plotter.axis.text(0.5,0.9,'Sum',transform=axis.transAxes)
