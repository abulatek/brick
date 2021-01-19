from spectral_cube import SpectralCube
import pyspeckit
# import pyspeckit.Spectrum as sp
from astroquery.splatalogue import Splatalogue

results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
cube_path = 'source_ab_138_spw31_clean_2sigma_n100000_masked_3sigma.image'

print("Now importing cube with spectral_cube")
scube = SpectralCube.read(results_path+cube_path, format='casa_image')
print("Successfully imported cube!")
print("Now transferring cube to pyspeckit")
pcube = pyspeckit.Cube(cube=scube, header=scube.hdulist)
print("Successfully transferred cube!")
print("Now getting spectrum from cube")
pspec = pcube.get_spectrum(314,199)
print("Successfully got spectrum!")

print("Now opening plotter")
pyspeckit.plotter()
print("Successfully opened plotter!")
print("Now fitting a spectrum")
pyspeckit.specfit(fittype='gaussian')
print("Successfully fit spectrum!")

# Splatalogue.query_lines(min_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 - dv/3e5)*u.Hz, max_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 + dv/3e5)*u.Hz, line_lists=['SLAIM']).pprint(max_lines=1000)
