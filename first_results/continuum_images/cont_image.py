# sbatch --job-name=cont_image_incl_adjparam --output=cont_image_incl_adjparam-%j.log --account=astronomy-dept --qos=astronomy-dept-b --ntasks=8 --mail-type=ALL --mail-user=abulatek@ufl.edu --nodes=1 --mem=24gb --time=96:00:00 --constraint=el8 --wrap "/blue/adamginsburg/adamginsburg/casa/casa-6.6.0-2-py3.8.el8/bin/casa -c 'execfile(\"/blue/adamginsburg/abulatek/brick/first_results/continuum_images/cont_image.py\")'"

## REMEMBER: you need to escape the quotes!! 

import numpy as np
import string
import os

# https://github.com/ALMA-IMF/reduction/blob/master/reduction/parse_contdotdat.py

def parse_contdotdat(filepath):
    selstr = ""
    selections = []
    with open(filepath, 'r') as fh:
        lines = fh.readlines()
    nlines = len(lines)
    spw = None
    for ii, line in enumerate(lines):
        if 'SpectralWindow' in line:
            spw = line.split()[1]
        if "LSRK" in line:
            selections.append(line.split()[0])
        if line.strip() == '' and spw is not None:
            selstr = selstr + f'{spw}:{";".join(selections)}'
            spw = None
            selections = []
            if ii < nlines - 1:
                selstr = selstr + ","

    return selstr

ms_prefixes = [
'/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a31/group.uid___A001_X1465_X3a32/member.uid___A001_X1465_X3a33', # Custom-made cont.dat file
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a41/group.uid___A001_X1465_X3a42/member.uid___A001_X1465_X3a43', # Custom-made cont.dat file
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a51/group.uid___A001_X1465_X3a52/member.uid___A001_X1465_X3a53',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a55/group.uid___A001_X1465_X3a56/member.uid___A001_X1465_X3a57',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a59/group.uid___A001_X1465_X3a5a/member.uid___A001_X1465_X3a5b',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a5d/group.uid___A001_X1465_X3a5e/member.uid___A001_X1465_X3a5f',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a61/group.uid___A001_X1465_X3a62/member.uid___A001_X1465_X3a63'
]

# window_names = ['data_31',
#                 'data_41',
#                 'data_51',
#                 'data_55',
#                 'data_59',
#                 'data_5d',
#                 'data_61',]


ms_dirs = [ms_prefix + '/calibrated/calibrated_final.ms' for ms_prefix in ms_prefixes]

cont_incls = [parse_contdotdat(ms_prefix + '/calibration/cont.dat') for ms_prefix in ms_prefixes]

# Mappings for colloquial names
mappings = {
    'data_31': 'X3a33',
    'data_41': 'X3a43',
    'data_51': 'X3a53',
    'data_55': 'X3a57',
    'data_59': 'X3a5b',
    'data_5d': 'X3a5f',
    'data_61': 'X3a63'
}

cellsize = {
    'data_31': '0.2arcsec',
    'data_41': '0.1arcsec',
    'data_51': '0.2arcsec',
    'data_55': '0.2arcsec',
    'data_59': '0.2arcsec',
    'data_5d': '0.2arcsec',
    'data_61': '0.2arcsec'
}

mappings_mirror = {value: key for key, value in mappings.items()}

image_prefix = '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_'
image_suffix = '_mtmfs_incl_adjparam'

image_paths = [image_prefix + mappings_mirror[ms_prefix.split('_')[-1]] + image_suffix for ms_prefix in ms_prefixes]

for (ms_dir, cont_incl, image_path) in zip(ms_dirs, cont_incls, image_paths):
    if not os.path.exists(image_path+'.image.tt0'):
        window_name = mappings_mirror[ms_dir.split('/')[-3].split("_")[-1]]
        print(ms_dir)
        print(cont_incl)
        print(image_path)
        cell = cellsize[window_name]
        print(cell)
        print(window_name)
        tclean(vis=ms_dir, imagename=image_path, specmode='mfs',
                deconvolver='mtmfs', nterms = 2, spw=cont_incl,
                field='BrickMaser', niter=50000, imsize=512, cell=cell,
                pbcor=True, usemask='pb', pbmask=0.18, cyclefactor=3, )
