from spectral_cube import SpectralCube, OneDSpectrum
import pyspeckit

file_dir = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/'
image_file = file_dir + 'BrickMaser_87_spw25.image.pbcor.fits'

cube = SpectralCube.read(image_file)
print(cube)

spectrum = cube[:,440,440]

sp = pyspeckit.Spectrum(data = spectrum) # May need to include data = spectrum, need xarr to be spectral axis

# The following will define the spectral axis units using the specified velocity convention
spectral_axis = sp_mJyperbeam.with_spectral_unit(u.GHz, velocity_convention=velconvention, rest_value=restf).spectral_axis
sp = pyspeckit.Spectrum(xarr=spectral_axis, data=cubespec_mJyperbeam.value, header={}) # Extract spectrum from cubespec_mJyperbeam


sp.plotter()
sp.plotter.savefig('basic_plot_example.png')

