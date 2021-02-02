from spectral_cube import SpectralCube
from astropy import units as u

file_dir = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/'
image_file = 'source_ab_138_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18' # no file ending

cube = SpectralCube.read(file_dir+image_file+'.image', format = 'casa_image')
cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=153.86508*u.GHz)
cube.write(file_dir+image_file+'_vel_HNCO.fits', format='fits', overwrite=True)
