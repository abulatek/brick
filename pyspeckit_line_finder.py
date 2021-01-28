# Manual single-component line fitting with pyspeckit

from spectral_cube import SpectralCube
import pyspeckit
from pyspeckit import units
from astroquery.splatalogue import Splatalogue

results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
cube_path = 'source_ab_138_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'

print("Now importing cube with spectral_cube")
cube = SpectralCube.read(results_path+cube_path, format='casa_image')
print("Shape of cube:",cube.shape)
print("Successfully imported cube!") 

'''
print("Now testing pyspeckit functionality")
print("Current pyspeckit implementation:")
spec = cube[:,200,300]
xarr_cur = units.SpectroscopicAxis(spec.spectral_axis)
print(type(xarr_cur))
print("Modified implementation:")
xarr_mod = cube.spectral_axis
print(type(xarr_mod))
'''

print("Now extracting spectrum with pyspeckit")
sp = pyspeckit.Spectrum(data=cube[:,200,300], xarr=cube.spectral_axis) # y, x (column, row)
print("Successfully transferred spectrum!")
print("Now opening plotter")
sp.plotter()
print("Successfully opened plotter!")
'''
print("Now fitting a spectrum")
pyspeckit.specfit(fittype='gaussian')
print("Successfully fit spectrum!")
'''
# Splatalogue.query_lines(min_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 - dv/3e5)*u.Hz, max_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 + dv/3e5)*u.Hz, line_lists=['SLAIM']).pprint(max_lines=1000)
