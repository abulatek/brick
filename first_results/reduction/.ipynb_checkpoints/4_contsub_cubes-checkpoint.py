# For parallelization, which helps to make convolution faster
import dask
dask.config.set(scheduler = 'threads', num_workers = 8) # Used to be 8 but I don't understand parallelization
# from dask.diagnostics import ProgressBar
# ProgressBar().register()
import warnings
warnings.filterwarnings('ignore')
import os

from astropy import units as u
from spectral_cube import SpectralCube

results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/smoothed_cubes_K_w_pbcor/'
destination = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/contsub_smoothed_cubes_K_w_pbcor/'

# Get names of spectral windows
# I think these are in vaguely the right order? Some might be repeats
freq_spws = ['87_spw25', '89_spw27', '91_spw25', '93_spw27', '95_spw25', '97_spw27', '98_spw29', '99_spw31', '102_spw23', '102_spw29', '104_spw25', '103_spw31', '106_spw29', '107_spw31', '110_spw29', '111_spw31', '112_spw27', '114_spw29', '127_spw65', '129_spw67', '130_spw105', '132_spw107', '134_spw45', '135_spw47', '137_spw85', '137_spw69', '139_spw71', '140_spw109', '142_spw111', '142_spw27', '144_spw49', '146_spw51', '147_spw89', '149_spw91', '142_spw27', '151_spw29', '152_spw31', '244_spw65', '245_spw67', '247_spw105', '249_spw107', '250_spw25', '252_spw27', '254_spw85', '255_spw87', '257_spw45', '259_spw47', '259_spw71', '261_spw109', '263_spw111', '264_spw29', '266_spw31', '268_spw89', '270_spw91', '271_spw49', '273_spw51']

for freq_spw in freq_spws:
    fn_orig = 'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pbcor.image.commonbeam.K.fits'
    path_orig = results+fn_orig
    outfile = destination+fn_orig.replace('.commonbeam.K.fits','.commonbeam.contsub.K.fits')
    if not os.path.exists(outfile):
        # Import cube and convert it to a common beam
        smoothed_cube = SpectralCube.read(path_orig, use_dask=True)
        smoothed_cube.allow_huge_operations=True
        smoothed_contsub_cube = smoothed_cube - smoothed_cube.median(axis=0) # Do continuum subtraction
        # Write the new cube
        smoothed_contsub_cube.write(outfile, format='fits', overwrite=True)
        print(f"Finished {fn_orig}")