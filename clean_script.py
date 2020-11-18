data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'
results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

output_name = results_path+'source_ab_138_spw25_clean_2sigma_masked_2sigma' # Do not add .image
mask_name = results_path+'source_ab_138_spw25_dirty_512_2sigma_e3_d4.mask'
threshold_val = '3.5mJy'
iterations = 10000

tclean(vis=data_path, imagename=output_name, field='2', spw='25', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold=threshold_val, mask=mask_name, niter=iterations, interactive=False, chanchunks=-1)

