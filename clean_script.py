data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'
results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

tclean(vis=data_path, imagename=results_path+'source_ab_138_spw25_clean_2sigma_masked', field='2', spw='25', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='3.5mJy', mask=results_path+'source_ab_138_spw25_dirty_512.mask', niter=10000, interactive=False, chanchunks=-1)

