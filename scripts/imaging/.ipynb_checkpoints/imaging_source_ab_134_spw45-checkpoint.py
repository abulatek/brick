# Script for imaging, meant to be run in CASA6 (necessary to access spectral_cube)
# Includes section for creating dirty map, masking dirty map for imaging (visualization elements removed), and imaging

# Import statements
import numpy as np
from astropy.stats import mad_std
from spectral_cube import SpectralCube
from astropy import units as u
from casatools import image
import scipy.ndimage
import os

# Define variables
data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'
output_dir = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
image_prefix = output_dir+'source_ab_134_spw45'
spw = '45'

def make_dirty_map(data_path, image_prefix, spw):
    '''Make a dirty map from a spectral window in a dataset. Parameters are pre-tuned but can be edited.'''
    dirty_map_image = image_prefix+'_dirty'
    print("Now creating dirty map:",dirty_map_image)
    tclean(vis=data_path, imagename=dirty_map_image, field='2', spw=spw, specmode='cube', gridder='standard',
	   cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='0.0mJy', interactive=False, chanchunks=-1)
    print("Dirty map creation successful!")
    return dirty_map_image # Return image name for use later

def measure_sigma(cube, nchan=50):
    '''Calculate a more robust version of the standard deviation in a cube with the median absolute deviation.
    Uses the middle 50 frequency channels by default, but this can be changed.'''
    print("Now measuring the mad_std in the dirty map")
    width = nchan/2
    freq_lo, freq_hi = cube.shape[0]/2-width, cube.shape[0]/2+width
    sigma = np.nanmedian(cube[int(freq_lo):int(freq_hi),:,:].mad_std(axis=0)) # Ignores NaN values, will be in mJy/beam
    sigma_mJy_beam = sigma.to(u.mJy/u.beam)
    print("Standard deviation calculation successful! Result:",round(sigma_mJy_beam.value,1),"mJy/beam")
    return round(sigma_mJy_beam.value,1)

def make_mask(dirty_cube, dirty_map_image, sigma, erosion_iter, dilation_iter, erosion_dilation=True):
    '''Masks a dirty cube at three times a provided sigma value.'''
    mask_image = dirty_map_image+'_3sigma' # Initialize mask filename
    print("Now masking the cube at 3 sigma")
    mask = dirty_cube > 3.*sigma*u.mJy/u.beam
    print("Masking successful!")
    print("Now computing boolmask")
    ia = image()
    boolmask = mask.include().compute()
    print("Boolmask computation successful!")
    print("Now getting coordinates from original image")
    ia.open(dirty_map_image+'.image')
    cs = ia.coordsys()
    ia.close()
    print("Coordinates retrieved successfully!")
    print("Now outputting the mask file")
    ia.fromarray(outfile=mask_image+'.mask', pixels=boolmask.astype('float')[:,None,:,:].T, csys=cs.torecord(), overwrite=True)
    ia.close()
    print("Mask file output successful!")
    if erosion_dilation == False:
        return mask_image
    if erosion_dilation == True:
        print("Now calculating the eroded/dilated boolmask")
        boolmask_e_d = scipy.ndimage.binary_dilation(scipy.ndimage.binary_erosion(boolmask, iterations=erosion_iter), iterations=dilation_iter)
        print("Eroded/dilated boolmask calculation successful!")
        print("Now defining the eroded/dilated mask to export")
        ia.open(mask_image+'.mask') # Open the mask we just created to be rewritten
        boolmask_e_d_export = boolmask_e_d[:,None,:,:].T.astype('int')
        print("Now outputting the eroded/dilated mask")
        ia.putchunk(pixels=boolmask_e_d_export)
        ia.close()
        print("Eroded/dilated mask output successful!")
        print("Now renaming mask file for consistency")
        mask_e_d_image = mask_image+'_e'+str(erosion_iter)+'_d'+str(dilation_iter)
        os.rename(mask_image+'.mask',mask_e_d_image+'.mask')
        print("Eroded/dilated mask file successfully renamed!")
        return mask_e_d_image

def deep_clean(data_path, output_dir, image_prefix, spw, sigma, nsigma, mask_image, niter):
    '''Perform a deep clean on a dataset.'''
    cleaned_image = image_prefix+'_clean_'+str(int(nsigma))+'sigma_n'+str(niter)+'_masked_3sigma_pbmask0p18' # Do not add .image because CASA just wants prefix
    print("Now beginning deep clean:",cleaned_image)
    tclean(vis=data_path, imagename=cleaned_image, field='2', spw=spw, specmode='cube', gridder='standard', 
           cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold=str(round(nsigma*sigma,1))+'mJy', mask=mask_image+'.mask',
           pbmask=0.18, niter=niter, interactive=False, chanchunks=-1)
    print("Deep clean complete!")
    return cleaned_image

dirty_map_image = make_dirty_map(data_path, image_prefix, spw)
dirty_cube = SpectralCube.read(dirty_map_image+'.image') # Not sure if this is best method, but helps us avoid multiple read-ins
sigma = measure_sigma(dirty_cube) # Will have units of mJy/beam, but be a float rounded to 1 decimal place
mask_image = make_mask(dirty_cube, dirty_map_image, sigma, 2, 2)
cleaned_image = deep_clean(data_path, output_dir, image_prefix, spw, sigma, 2., mask_image, 50000)
