## Script for generating primary beam-corrected cubes

## Import statements
import glob
from spectral_cube import SpectralCube
import astropy.units as u
import sys
import os
# For parallelization, which helps to make convolution faster
import dask
dask.config.set(scheduler = 'threads', num_workers = 8)
from dask.diagnostics import ProgressBar
    
## Grab spectral window from command line arguments
freq_spw = sys.argv[1] # Will be a string; e.g. 87_spw25
if len(sys.argv) > 2:
    progressbar = sys.argv[2] == "True"
else:
    progressbar = False
    
if progressbar:
    ProgressBar().register()

## Get name of cube and the pb file for that cube
pb_fn = glob.glob(f"/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab*{freq_spw}*clean*.pb")
cube_fn = glob.glob(f"/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab*{freq_spw}*clean*pbmask0p18.image") # We need the units of this... and no one knows why

## Generate output filename
outfile = cube_fn[0].replace(".image",".pbcor.image.fits") # THIS SHOULD BE .pbcor.image.fits

if not os.path.exists(outfile):
    pb = SpectralCube.read(pb_fn[0], format = 'casa_image')
    cube = SpectralCube.read(cube_fn[0], format = 'casa_image')

    ## Use division to correct the cube, and deal with the units issue
    cube_pbcor = cube.unitless/(pb.unitless)
    cube_pbcor_unit = cube_pbcor * (u.Jy/u.beam)

    ## Write the pb-corrected cube
    print(outfile)
    cube_pbcor_unit.write(outfile, format='fits')
    # AFTER THIS, I WENT AND MANUALLY RENAMED ALL OF THESE OUTPUT FILES TO END IN .pbcor.image.fits
    # BECAUSE THE STRING OUTFILE DOESN'T CONTAIN A .fits AT THE END. BUT THEN I CLEANED A FEW MORE
    # CASA IMAGES, AND USED THIS FILE FOR THEM
