# Written by Adam Ginsburg (from his Github page)

import glob
import pyspeckit
import pylab as pl
from spectral_cube import OneDSpectrum
from astropy.io import fits
from astropy import units as u

pl.close(1)
pl.close(2)
fig=pl.figure(1, figsize=(25,10))

spectra = pyspeckit.Spectra(pyspeckit.Spectrum(x) for x in glob.glob("spectra/*max.fits"))
spectra.plotter(figure=fig)

spectra.plotter.figure.savefig('full_spectrum_max.png', dpi=200, bbox_inches='tight')



fig2=pl.figure(2, figsize=(25,10))

kspectra = [OneDSpectrum.from_hdu(fits.open(x)).to(u.K) for x in glob.glob("spectra/*max.fits")]
kspectra_ps = [pyspeckit.Spectrum.from_hdu(kspec.hdu) for kspec in kspectra]
spectra_K = pyspeckit.Spectra(kspectra_ps)
spectra_K.plotter(figure=fig2)

spectra_K.plotter.figure.savefig('full_spectrum_max_K.png', dpi=200, bbox_inches='tight')
spectra_K.plotter.figure.savefig('full_spectrum_max_K.pdf', dpi=200, bbox_inches='tight')
