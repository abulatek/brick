data_path = '/blue/adamginsburg/abulatek/brick/symlinks/data_31/calibrated_final.ms'

tclean(vis=data_path, imagename='/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_152_spw31_dirty_512', field='2', spw='31', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='0.0mJy', interactive=False, chanchunks=-1)
