# sbatch --job-name=cont_image_incl_adjparam_split --output=cont_image_incl_adjparam_split-%j.log --account=astronomy-dept --qos=astronomy-dept-b --ntasks=8 --mail-type=ALL --mail-user=abulatek@ufl.edu --nodes=1 --mem=24gb --time=96:00:00 --constraint=el8 --wrap "/blue/adamginsburg/adamginsburg/casa/casa-6.6.0-2-py3.8.el8/bin/casa -c 'execfile(\"/blue/adamginsburg/abulatek/brick/first_results/continuum_images/cont_image_splitB4B6.py\")'"

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
'/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a41/group.uid___A001_X1465_X3a42/member.uid___A001_X1465_X3a43', # Custom-made cont.dat file
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a51/group.uid___A001_X1465_X3a52/member.uid___A001_X1465_X3a53',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a55/group.uid___A001_X1465_X3a56/member.uid___A001_X1465_X3a57',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a59/group.uid___A001_X1465_X3a5a/member.uid___A001_X1465_X3a5b',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a5d/group.uid___A001_X1465_X3a5e/member.uid___A001_X1465_X3a5f',
# '/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a61/group.uid___A001_X1465_X3a62/member.uid___A001_X1465_X3a63'
]

ms_dirs = [ms_prefix + '/calibrated/calibrated_final.ms' for ms_prefix in ms_prefixes]

# cont_incl = [parse_contdotdat(ms_prefix + '/calibration/cont.dat') for ms_prefix in ms_prefixes] # Want to do DIY cont.dats so I can sep
cont_incl = ['65:125.083744140625~126.958255859375GHz,67:126.729244140625~127.330000000000GHz;126.390000000000~127.820000000000GHz;127.880000000000~128.603755859375GHz,105:128.46899414052498~128.73000000000000GHz;128.80000000000000~129.10000000000000GHz;129.15000000000000~129.33000000000000GHz;129.38000000000000~130.22000000000000GHz;130.29000000000000~130.34350585927498GHz,107:130.17211914052498~130.22000000000000GHz;130.29000000000000~131.84000000000000GHz;131.92000000000000~132.04663085927498GHz,45:131.92000000000000~132.85000000000000GHz;132.92000000000000~133.74975585927498GHz,47:133.57836914052498~133.75000000000000GHz;133.80000000000000~135.27000000000000GHz;135.31000000000000~135.45288085927498GHz,85:135.30000000000000~136.420000000000GHz;136.48000000000000~136.700000000000GHz;136.75000000000000~137.15600585927498GHz,87:136.98461914052498~137.15600585927498GHz;138.729244140625~138.85913085927498GHz,69:137.146244140625~137.330000000000GHz;137.390000000000~138.130000000000GHz;138.220000000000~138.710000000000GHz;138.730000000000~139.020755859375GHz,71:138.729244140625~139.45000000000000GHz;139.50000000000000~140.603755859375GHz,25:138.68774414052498~139.45000000000000GHz;139.50000000000000~140.56225585927498GHz','27:140.39086914052498~140.80000000000000GHz;140.87000000000000~141.54000000000000GHz;141.62000000000000~142.26538085927498GHz,109:140.46899414072502~140.80000000000000GHz;140.87000000000000~141.54000000000000GHz;141.62000000000000~142.34350585947502GHz,111:142.17211914072502~143.81000000000000GHz;143.89000000000000~144.04663085947502GHz,49:143.87524414072502~144.58000000000000GHz;144.63000000000000~144.82000000000000GHz;144.84000000000000~145.05000000000000GHz;145.15000000000000~145.52000000000000GHz;145.62000000000000~145.74975585947502GHz,51:145.65000000000000~145.90000000000000GHz;145.95000000000000~146.33000000000000GHz;146.39000000000000~146.60000000000000GHz;146.64000000000000~146.89000000000000GHz;147.00000000000000~147.05000000000000GHz;147.20000000000000~147.45288085947502GHz,89:147.28149414072502~149.15600585947502GHz,91:148.98461914072502~150.45000000000000GHz;150.53000000000000~150.78000000000000GHz,29:150.68774414072502~150.80000000000000GHz;150.87000000000GHz~152.56225585947502GHz,31:152.39086914072502~152.68000000000000GHz;152.71000000000000~153.76000000000000GHz;153.88000000000000~154.26538085947502GHz','65:241.900000000000~243.150000000000GHz;243.200000000000~243.561630859275GHz,67:243.390244140525~243.850000000000GHz;243.900000000000~244.000000000000GHz;244.020000000000~244.200000000000GHz;244.230000000000~244.640000000000GHz;244.720000000000~244.850000000000GHz;244.950000000000~245.080000000000GHz;245.120000000000~245.264755859275GHz,105:245.093369140525~245.550000000000GHz;245.600000000000~246.967880859275GHz,107:246.796494140525~247.200000000000GHz;247.220000000000~248.150000000000GHz;248.230000000000~248.671005859275GHz,25:248.499619140525~250.374130859275GHz,27:250.202744140525~251.750000000000GHz;251.900000000000~252.077255859275GHz,85:251.90586914052500~251.94000000000000GHz;251.96000000000000~252.03000000000000GHz;252.07000000000000~252.18000000000000GHz;252.22000000000000~253.50000000000000GHz;253.56000000000000~253.780380859275GHz,87:253.608994140525~253.800000000000GHz;254.200000000000~254.620000000000GHz;254.720000000000~255.300000000000GHz;255.400000000000~255.483505859275GHz,45:255.312119140525~257.186630859275GHz,47:257.015244140525~257.420000000000GHz;257.500000000000~258.889755859275GHz,69:257.71296529999995~258.0924130000GHz;258.15243009999995~258.1912450000GHz;258.25124500000000~258.93129480000005GHz','69:259.00521690000000~259.56163085947503GHz,71:259.39024414072503~260.20000000000000GHz;260.23000000000000~260.45000000000000GHz;260.55000000000000~261.26475585947503GHz,109:261.09336914072503~261.2200000000GHz;261.25000000000000~261.7500000000GHz;261.82000000000000~261.9000000000GHz;262.1000000000~262.96788085947503GHz,111:262.79649414072503~263.7200000000GHz;263.7800000000~264.67100585947503GHz,29:264.49961914072503~265.2000000000GHz;265.30000000000000~265.6500000000GHz;266.00000000000000~266.1500000000GHz;266.2500000000~266.37413085947503GHz,31:266.20274414072503~266.7500000000GHz;266.85000000000000~267.5000000000GHz;267.5800000000~268.07725585947503GHz,89:267.90586914072503~269.78038085947503GHz,91:269.60899414072503~270.4500000000GHz;270.5000000000~271.48350585947503GHz,49:271.31211914072503~272.8200000000GHz;272.0000000000~273.18663085947503GHz,51:273.01524414072503~274.4000000000GHz;274.500000000000~274.889755859475GHz']

# Mappings for colloquial names
mappings = {
    'data_31': 'X3a33',
    'data_41': 'X3a43',
    # 'data_51': 'X3a53',
    # 'data_55': 'X3a57',
    # 'data_59': 'X3a5b',
    # 'data_5d': 'X3a5f',
    # 'data_61': 'X3a63'
}

cellsize = {
    'data_31': '0.2arcsec',
    'data_41': '0.1arcsec',
    # 'data_51': '0.2arcsec',
    # 'data_55': '0.2arcsec',
    # 'data_59': '0.2arcsec',
    # 'data_5d': '0.2arcsec',
    # 'data_61': '0.2arcsec'
}

mappings_mirror = {value: key for key, value in mappings.items()}

image_prefix = '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_'
image_suffix = '_mtmfs_incl_adjparam_split'

image_paths = [image_prefix + mappings_mirror[ms_prefix.split('_')[-1]] + image_suffix for ms_prefix in ms_prefixes]

i=0
for (ms_dir, image_path) in zip(ms_dirs, image_paths):
    if not os.path.exists(image_path+'.image.tt0'):
        window_name = mappings_mirror[ms_dir.split('/')[-3].split("_")[-1]]
        print(ms_dir)
        print(cont_incl)
        print(image_path)
        cell = cellsize[window_name]
        print(cell)
        print(window_name)
        tclean(vis=ms_dir, imagename=image_path+'_first', specmode='mfs',
                deconvolver='mtmfs', nterms = 2, spw=cont_incl[i],
                field='BrickMaser', niter=50000, imsize=512, cell=cell,
                pbcor=True, usemask='pb', pbmask=0.18, cyclefactor=3, )
    if not os.path.exists(image_path+'.image.tt0'):
        window_name = mappings_mirror[ms_dir.split('/')[-3].split("_")[-1]]
        print(ms_dir)
        print(cont_incl)
        print(image_path)
        cell = cellsize[window_name]
        print(cell)
        print(window_name)
        tclean(vis=ms_dir, imagename=image_path+'_second', specmode='mfs',
                deconvolver='mtmfs', nterms = 2, spw=cont_incl[i+1],
                field='BrickMaser', niter=50000, imsize=512, cell=cell,
                pbcor=True, usemask='pb', pbmask=0.18, cyclefactor=3, )
    i+=2
