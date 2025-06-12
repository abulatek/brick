data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'
results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

output_name = results_path+'source_ab_151_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18' # Do not add .image
mask_name = results_path+'source_ab_151_spw29_dirty_512_3sigma_e2_d2.mask'
threshold_val = '4.24mJy' # needs to be changed
iterations = 50000

tclean(vis=data_path, imagename=output_name, field='2', spw='29', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold=threshold_val, mask=mask_name, pbmask=0.18, niter=iterations, interactive=False, chanchunks=-1)
