from metadata_tools import determine_phasecenter, determine_imsize

vis = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a31/group.uid___A001_X1465_X3a32/member.uid___A001_X1465_X3a33/calibrated/calibrated_final.ms'

determine_phasecenter(ms=vis,field=2)
#determine_imsize(ms=vis)
