full_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'

# tclean(vis=full_path, imagename='/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_138_spw25_dirty_512', field='2', spw='25', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='0.0mJy', interactive=False, chanchunks=-1)

tclean(vis=full_path, imagename='source_ab_138_spw25_clean_2sigma_masked', field='2', spw='25', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='3.5mJy', mask='source_ab_138_spw25_dirty_512.mask', niter=10000, interactive=False, chanchunks=-1)

