from scousepy import scouse
from astropy.io import fits
import os
import matplotlib.pyplot as pl
pl.ion()

def run_scousepy():
  datadirectory =  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
  outputdir     =  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/scouse_output/simple_example_run/'
  filename      =  'source_ab_138_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18_vel_HNCO'
  ppv_vol       = [0.0,50.0,170,240,270,340]
  wsaa          = [8.0]
  tol           = [3.0, 1.0, 3.0, 3.0, 0.5]
  verb          = True
  fittype       = 'gaussian'
  njobs         = 4
  mask          = 0.3

  #==========================================================================#
  # Stage 1
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_1/s1.scousepy'):
      s = scouse(outputdir=outputdir, filename=filename, fittype=fittype,
                 datadirectory=datadirectory)
      s.load_stage_1(outputdir+filename+'/stage_1/s1.scousepy')
      s.load_cube(fitsfile=datadirectory+filename+".fits")
  else:
      s = scouse.stage_1(filename, datadirectory, wsaa,
                         ppv_vol=ppv_vol,
                         outputdir=outputdir,
                         mask_below=mask,
                         fittype=fittype,
                         verbose = verb,
                         write_moments=True,
                         save_fig=True)

  #==========================================================================#
  # Stage 2
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_2/s2.scousepy'):
      s.load_stage_2(outputdir+filename+'/stage_2/s2.scousepy')
  else:
      s = scouse.stage_2(s, verbose=verb, write_ascii=True)

  #==========================================================================#
  # Stage 3
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_3/s3.scousepy'):
      s.load_stage_3(outputdir+filename+'/stage_3/s3.scousepy')
  else:
      s = scouse.stage_3(s, tol, njobs=njobs, verbose=verb)

  #==========================================================================#
  # Stage 4
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_4/s4.scousepy'):
      s.load_stage_4(outputdir+filename+'/stage_4/s4.scousepy')
  else:
      s = scouse.stage_4(s, verbose=verb)

  #==========================================================================#
  # Stage 5
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_5/s5.scousepy'):
      s.load_stage_5(outputdir+filename+'/stage_5/s5.scousepy')
  else:
      s = scouse.stage_5(s, blocksize = 6,
                            figsize = [18,10],
                            plot_residuals=True,
                            verbose=verb)

  #==========================================================================#
  # Stage 6
  #==========================================================================#
  if os.path.exists(outputdir+filename+'/stage_6/s6.scousepy'):
      s.load_stage_6(outputdir+filename+'/stage_6/s6.scousepy')
  else:
      s = scouse.stage_6(s, plot_neighbours=True,
                            radius_pix = 2,
                            figsize = [18,10],
                            plot_residuals=True,
                            write_ascii=True,
                            verbose=verb)

  return s

s = run_scousepy()
