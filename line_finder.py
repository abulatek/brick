from spectral_cube import SpectralCube
import pyspeckit as sp
from astroquery.splatalogue import Splatalogue

results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
cube_path = 'source_ab_138_spw31_clean_2sigma_n100000_masked_3sigma.image'

cube = SpectralCube.read(results_path+cube_path, format='casa_image')
pcube = pyspeckit.Cube(cube)
pspec = pcube.get_spectrum(314,199)

sp.plotter()
sp.specfit(fittype='gaussian')

# Splatalogue.query_lines(min_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 - dv/3e5)*u.Hz, max_frequency=sp.specfit.parinfo['SHIFT0'] * (1 - velocity/3e5 + dv/3e5)*u.Hz, line_lists=['SLAIM']).pprint(max_lines=1000)
