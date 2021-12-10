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
# cube = SpectralCube.read('temp_cube.fits').spectral_slab(218.14*u.GHz, 218.54*u.GHz) # orig version


# freq_lo = 146.8*u.GHz # 110.02*u.GHz # 146.8*u.GHz
# freq_hi = 147.2*u.GHz # 111.26*u.GHz# 147.2*u.GHz

# # alyssa is trying something again
# cube = SpectralCube.read('../methyl_cyanide/ch3cn_0_masked.fits')
# cube = SpectralCube.read('/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/BrickMaser_147_spw89.image.pbcor.fits').spectral_slab(freq_lo, freq_hi)

spw = '110_spw29' # '146_spw51'
cube_filename = f'/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/source_ab_{spw}_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(cube_filename, format='casa_image') #.spectral_slab(freq_lo, freq_hi) # was using 146_spw51. 110_spw29

print(cube)
# cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=147.1745883*u.GHz)
# cube = cube.with_spectral_unit(u.GHz, velocity_convention='radio', rest_value=147.1745883*u.GHz) # , rest_value=147.1745883*u.GHz
freq_axis = cube.spectral_axis

# this isn't really right b/c of beam issues, but it's Good Enough
# meanspec = cube[:,256,256]
# meanspec = cube[:,203,303] # max temp?


# yc, xc = 196,291
# meanspec = cube[:,yc,xc] # This is not really the meanspec
meanspec = cube.max(axis=(1,2)) #, how='slice', progressbar=True)
print(meanspec)


#import temp map, col density map, and mom1 map for inputting parameters
from astropy.io import fits
prefix = '/blue/adamginsburg/abulatek/brick/first_results/'
temp = fits.getdata(f'{prefix}temperature_map_test.fits')
log_N = fits.getdata(f'{prefix}log_N_test.fits')

masked_cube = SpectralCube.read(f'{prefix}methyl_cyanide/ch3cn_0_masked.fits', format='fits')
mom1 = masked_cube.moment1()

print(temp.shape, log_N.shape, mom1.shape)


# have to hack this, which is _awful_ and we need to fix it
meanspec_K = meanspec.value * cube.jtok_factors()

import pyspeckit

sp = pyspeckit.Spectrum(xarr=cube.spectral_axis, data=meanspec_K)


# alyssa problem solving?
# freq_lo = 146.8*u.GHz
# freq_hi = 147.2*u.GHz
# freq_lo = 100*u.GHz
# freq_hi = 105*u.GHz
# xarr = pyspeckit.spectrum.units.SpectroscopicAxis(np.linspace(freq_lo,freq_hi,100), units='GHz')
# xarr = np.linspace(freq_lo,freq_hi,1000) # this is bad no use


# # generate a random spectrum
# sigma = 10.*u.GHz
# center = (freq_lo + freq_hi)/2.
# synth_data = np.exp(-(xarr-center)**2/(sigma**2 * 2.))*u.GHz

# # Add noise
# stddev = 0.1*u.GHz
# noise = np.random.randn(xarr.size)*stddev
# error = stddev*(np.ones_like(synth_data).value)
# data = noise+synth_data

# print('xarr:',xarr.unit)
# print('sigma:',sigma.unit)
# print('center:',center.unit)
# print('synth_data:',synth_data.unit)
# print('stddev:',stddev.unit)
# print('noise:',noise.unit)
# print('error:',error.unit)
# print('data:',data.unit)

# sp = pyspeckit.Spectrum(data=data, error=error, xarr=xarr,
#                         xarrkwargs={'unit':'GHz'},
#                         unit='K')

# # alyssa is testing
# results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/spectra/'
# # freq_spw = '139_spw71'
# freq_spw = '146_spw51'
# max_fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
# # mean_fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
# from spectral_cube import OneDSpectrum
# from astropy.io import fits
# kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn)).to(u.K)
# sp = pyspeckit.Spectrum.from_hdu(kspectrum.hdu)
# sp.xarr.convert_to_unit('GHz')


# The LTE model doesn't include an explicit filling factor, so we impose one
fillingfactor = 0.05

# there's continuum in our spectrum, but it's also not modeled
offset = np.nanpercentile(sp.data.filled(np.nan), 20)

# look at a couple species
# Be cautious! Not all line catalogs have all of these species, some species
# can result in multiple matches from splatalogue, and definitely not all
# species have the same column and excitation
# species_list = ('CH3OH', 'HCCCN', 'H2CO', 'CH3OCHO', 'CH3OCH3', 'C2H5CN',) # 'CN, v = 0, 1', 'H2CS',) #CNv=0?
# species_list = ('H2CS',)
species_list = ('HNCO', 'CH3CN',) # 'C18O', '13CO', #, 'HNCO'
# REMEMBER THE COMMA !

# molecules that didn't work: 'H13CN', 'H13CO+', 'HN13C', 'C18O', '13CO',
# I think I don't know how to search for these, and I don't know which database I should be picking names from. Need to find that out

# make a set of subplots:
import pylab as pl
# pl.close(1)
fig, axes = pl.subplots(len(species_list)+1, 1, figsize=(18,14), num=1)
pl.draw()

from pyspeckit.spectrum.models import lte_molecule
mods = []
for species, axis in zip(species_list, axes):
    print(species, axis)
    freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters_JPL(species, fmin=sp.xarr.min(), fmax=sp.xarr.max()) # this was the original version; it works with a working spectrum
    
    
    # now we need to improvise since alyssa didn't get a methyl cyanide cube oops
#     freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters_JPL(species, fmin=freq_lo, fmax=freq_hi) # alyssa's edited version
    
    
    #freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(species, fmin=200*u.GHz, fmax=250*u.GHz, export_limit=1e5, molecule_name_jpl=species_jpl,
    #                                                                      line_lists=['SLAIM'])
#     mod = lte_molecule.generate_model(sp.xarr, mom1[yc, xc]*u.km/u.s, 1.5*u.km/u.s, temp[yc, xc]*u.K, 10**(log_N[yc, xc])*u.cm**-2, freqs, aij, deg, EU, partfunc) # original version
    mod = lte_molecule.generate_model(sp.xarr, 50*u.km/u.s, 1.5*u.km/u.s, 100*u.K, 5e15*u.cm**-2, freqs, aij, deg, EU, partfunc) + lte_molecule.generate_model(sp.xarr, 20*u.km/u.s, 1.5*u.km/u.s, 100*u.K, 1e16*u.cm**-2, freqs, aij, deg, EU, partfunc)
#     mod = lte_molecule.generate_model(xarr, 50*u.km/u.s, 3*u.km/u.s, 100*u.K, 1e14*u.cm**-2, freqs, aij, deg, EU, partfunc) # alyssa's edited version
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
