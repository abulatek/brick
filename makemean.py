# Adapted from Adam's makemax.py

#import os
#import glob
from spectral_cube import SpectralCube, StokesSpectralCube
#from spectral_cube.dask_spectral_cube import DaskSpectralCube

#import time

input_dir = "/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/"
input_file = "BrickMaser_87_spw25"
suffix = ".image.pbcor.fits"

output_dir = input_dir + "spectra/center"
output_file = "BrickMaser_87_spw25_emission"

cube = StokesSpectralCube.read(input_dir + input_file + suffix)
print(cube)

subcube = 0 # cube defined by DS9 region
mnspec = subcube.mean(axis=(1,2))

mnspec.write(output_dir + output_file + suffix, overwrite = True)

# Below this line is the old code that did not output a file

#dask = True

#if dask:
#	from dask.diagnostics import ProgressBar
#	pbar = ProgressBar()
#	pbar.register()

#for suffix in ("image.pbcor.fits",):
	#for fn in glob.glob(f"BrickMaser*{suffix}"):
#	for fn in glob.glob(f"BrickMaser_87_spw25"): # Just testing on one cube at a time
#		t0 = time.time()
#		outf_ = fn.replace(suffix,"center.fits")
#		outfn = f'spectra/center/{outf_}'
#		print(fn, outfn)
#		if not os.path.exists(outfn):
#			if not dask:
#				cube = SpectralCube.read(fn)
#				print(cube)
#				subcube = cube[:, 380:420, 380:420]  
#				mnspec = subcube.mean(axis=(1,2))
				#mxspec = cube.max(axis=(1,2), how='slice', progressbar=True)
#			else:
#				print(fn)
#				cube = DaskSpectralCube.read(fn.replace(".fits", ""), format='casa_image')
#				print(cube)
#				subcube = cube[:, 380:420, 380:420]  
#				mnspec = subcube.mean(axis=(1,2))
				#mxspec = cube.max(axis=(1,2))

#			mnspec.write(outfn,
#					overwrite=True)
#			print(time.time() - t0)
