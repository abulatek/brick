# Clean commands for 107 and 108 GHz methanol lines

vis = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a59/group.uid___A001_X1465_X3a5a/member.uid___A001_X1465_X3a5b/calibrated/calibrated_final.ms/'

# 107 GHz line
tclean(vis=vis, imagename='brick_107GHz_dasar_line_1mJy_multiscale_cyclefactor4_amt_negativethreshold0p5', field='BrickMaser', spw='29', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='1mJy', pbmask=0.18, niter=1000000, interactive=False, start = '106.988171 GHz', nchan = 48, deconvolver = 'multiscale', scales = [0, 3, 9, 27], usemask = 'auto-multithresh', negativethreshold = 0.5, cyclefactor = 4)

# 108 GHz line
tclean(vis=vis, imagename='brick_108GHz_dasar_line_1mJy_multiscale_cyclefactor4', field='BrickMaser', spw='31', specmode='cube', gridder='standard', cell=['0.2arcsec'], imsize=[512,512], weighting='natural', threshold='1mJy', pbmask=0.18, niter=1000000, interactive=False, cyclefactor = 4, start = '108.862921 GHz', nchan = 48, deconvolver = 'multiscale', scales = [0, 3, 9, 27])