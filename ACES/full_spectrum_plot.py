import glob
import pyspeckit
import pylab as pl
from spectral_cube import OneDSpectrum
from astropy.io import fits
from astropy import units as u

def plot_spectrum(fn):
#     pl.close(1)
    pl.close(2)
    
#     fig=pl.figure(1, figsize=(10,5))
#     spectrum = pyspeckit.Spectrum(fn)
#     spectrum.plotter(figure=fig)
#     pl.show()

    fig2=pl.figure(2, figsize=(10,5))
    kspectrum = OneDSpectrum.from_hdu(fits.open(fn)).to(u.K)
    kspectrum_ps = pyspeckit.Spectrum.from_hdu(kspectrum.hdu)
    kspectrum_ps.plotter(figure=fig2)
    pl.show()

def plot_spectrum_line_ids(fn, ids):
#     pl.close(1)
    pl.close(2)
    
#     fig=pl.figure(1, figsize=(12,5))
#     spectrum = pyspeckit.Spectrum(fn)
#     spectrum.plotter(figure=fig)
#     spectrum.plotter.line_ids(figure=fig, line_names=ids['Species'], # ChemicalName 
#                               line_xvals=ids['Freq'], xval_units='GHz', plot_kwargs={'color':'silver'})
#     pl.show()
    
    fig2=pl.figure(2, figsize=(12,5))
    kspectrum = OneDSpectrum.from_hdu(fits.open(fn)).to(u.K)
    kspectrum_ps = pyspeckit.Spectrum.from_hdu(kspectrum.hdu)
    kspectrum_ps.plotter(figure=fig2)
    kspectrum_ps.plotter.line_ids(figure=fig2, line_names=ids['Species'], # ChemicalName 
                                  line_xvals=ids['Freq'], xval_units='GHz', plot_kwargs={'color':'silver'})
    pl.show()
    
def plot_spectrum_line_ids_final(fn, ids, save=False):    
#     pl.close(1)
    pl.close(2)
    
#     fig=pl.figure(1, figsize=(12,5))
#     spectrum = pyspeckit.Spectrum(fn)
#     spectrum.plotter(figure=fig)
#     spectrum.plotter.line_ids(figure=fig, line_names=ids['Species'], # ChemicalName 
#                               line_xvals=ids['Freq'], xval_units='GHz', plot_kwargs={'color':'silver'})
#     pl.show()
    
    fig2=pl.figure(2, figsize=(12,5))
    kspectrum = OneDSpectrum.from_hdu(fits.open(fn)).to(u.K)
    kspectrum_ps = pyspeckit.Spectrum.from_hdu(kspectrum.hdu)
    kspectrum_ps.xarr.convert_to_unit('GHz')
    kspectrum_ps.plotter(figure=fig2)
    kspectrum_ps.plotter.line_ids(figure=fig2, line_names=ids['Species'], # ChemicalName 
                                  line_xvals=ids['Frequency (GHz)'], xval_units='GHz', 
                                  plot_kwargs={'color':'silver'})
    pl.show()
    if save == True:
        filename = ids['freq_spw'][0]+'_'+ids['Spectrum type'][0]+'spectrum.png'
        kspectrum_ps.plotter.figure.savefig(filename, dpi=200, bbox_inches='tight')