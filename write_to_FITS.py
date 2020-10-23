from spectral_cube import SpectralCube

file_dir = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/'
image_file = file_dir + 'source_ab_138_spw25_dirty.image'

cube = SpectralCube.read(image_file, format = 'casa_image')

cube.write(image_file+'.fits', format='fits')
