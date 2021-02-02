# Script for imaging, meant to be run in CASA6 (necessary to access spectral_cube)
# Includes section for creating dirty map, masking dirty map for imaging, and imaging

# Import statements
import numpy as np
from astropy.stats import mad_std
from spectral_cube import SpectralCube
from astropy import units as u

# Define variables
data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'
output_dir = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
dirty_map_image = output_dir+'source_ab_141_spw25_dirty_512'
spw = '25'

def make_dirty_map(data_path, dirty_map_image, spw):
	'''Make a dirty map from a spectral window in a dataset. Parameters are pre-tuned but can be edited.'''
	tclean(vis=data_path, imagename=dirty_map_image, field='2', spw=spw, specmode='cube', gridder='standard',
	       cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='0.0mJy', interactive=False, chanchunks=-1)

def measure_sigma(dirty_map_image, nchan=50):
	'''Calculate a more robust version of the standard deviation in a cube with the median absolute deviation.
	Uses the middle 50 frequency channels by default, but this can be changed.'''
	cube = SpectralCube.read(dirty_map_image+'.image') # Want image file
	width = nchan/2 
	freq_lo, freq_hi = int(cube.shape[0]/2)-width, int(cube.shape[0]/2)+width
	sigma = np.nanmedian(cube[freq_lo:freq_hi,:,:].mad_std(axis=0)) # Ignores NaN values
	return sigma

def make_mask():
	return 1 # Work in progress

make_dirty_map(data_path, dirty_map_image, spw)
sigma = measure_sigma(dirty_map_image)

