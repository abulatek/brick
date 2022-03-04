########################################################
# Started Logging At: 2022-02-09 16:56:09
########################################################
########################################################
# # Started Logging At: 2022-02-09 16:56:10
########################################################
import time
start = time.time()

import pylab as pl
display_dpi = 150
# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# Get H2CS (template molecule) cube
freq_spw = '135_spw47'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
h2cscube = SpectralCube.read(fn, format='casa_image')
# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
# Get H2CS (template molecule) cube
freq_spw = '135_spw47'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
h2cscube = SpectralCube.read(fn, format='casa_image')
h2cscube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 133582251086.390 Hz:135456817280.560 Hz
from astroquery.splatalogue import Splatalogue
h2cscube.find_lines(chemical_name='H2CS', line_lists=['JPL'], 
                    show_upper_degeneracy=True).show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                show_upper_degeneracy=True, show_qn_code=True)
ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
ch3cntbl = ch3cntbl[::-1]
ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
ch3cntbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
h2cscube = h2cscube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                       rest_value=135.297811*u.GHz).spectral_slab(-10*u.km/u.s, 
                                                                                  80*u.km/u.s)
print(h2cscube)
h2cssubcube = h2cscube.spectral_slab(-10*u.km/u.s, 30*u.km/u.s)
def get_targetcube(comp):
    targetcube = ch3cncube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                             rest_value=ch3cn_freqs[comp]*u.GHz).spectral_slab(-10*u.km/u.s,
                                                                                               80*u.km/u.s)
    print(targetcube)
    targetsubcube = targetcube.spectral_slab(-10*u.km/u.s, 30*u.km/u.s)
    return targetcube, targetsubcube
targetcube, targetsubcube = get_targetcube(0)
h2cssubcube.max(axis=0).quicklook()
# This first quicklook has to be re-run, as it doesn't run the first time (not sure why)
targetsubcube.max(axis=0).quicklook()
# h2cssubcube = h2cssubcube.to(u.K)
# targetsubcube = targetsubcube.to(u.K)
h2cssubcube.beams, targetsubcube.beams
#[Out]# (<Beams [5.54782914e-11, 5.54786541e-11, 5.54790393e-11, 5.54794101e-11,
#[Out]#          5.54797749e-11, 5.54801437e-11, 5.54804993e-11, 5.54808641e-11,
#[Out]#          5.54812227e-11, 5.54815853e-11, 5.54819695e-11, 5.54823343e-11,
#[Out]#          5.54827113e-11, 5.54830658e-11, 5.54834461e-11, 5.54838098e-11,
#[Out]#          5.54841798e-11, 5.54845547e-11, 5.54849266e-11, 5.54852862e-11,
#[Out]#          5.54856530e-11, 5.54860146e-11, 5.54863998e-11, 5.54867451e-11,
#[Out]#          5.54871243e-11, 5.54874675e-11, 5.54878567e-11, 5.54882113e-11,
#[Out]#          5.54885853e-11, 5.54889633e-11, 5.54893137e-11, 5.54896959e-11,
#[Out]#          5.54900586e-11, 5.54904131e-11, 5.54907831e-11, 5.54911520e-11,
#[Out]#          5.54915250e-11, 5.54918918e-11] sr>,
#[Out]#  <Beams [4.77419483e-11, 4.77416592e-11, 4.77413732e-11, 4.77410841e-11,
#[Out]#          4.77407998e-11, 4.77405117e-11, 4.77402069e-11, 4.77399284e-11,
#[Out]#          4.77396318e-11, 4.77393506e-11, 4.77390635e-11, 4.77387792e-11,
#[Out]#          4.77384960e-11, 4.77381966e-11, 4.77379229e-11, 4.77376342e-11,
#[Out]#          4.77373499e-11, 4.77370526e-11, 4.77367591e-11, 4.77364635e-11,
#[Out]#          4.77361545e-11, 4.77358760e-11, 4.77355805e-11, 4.77352838e-11,
#[Out]#          4.77350102e-11, 4.77347108e-11, 4.77344296e-11, 4.77341388e-11,
#[Out]#          4.77338566e-11, 4.77335819e-11, 4.77332911e-11, 4.77330127e-11,
#[Out]#          4.77327332e-11, 4.77324537e-11, 4.77321619e-11, 4.77318701e-11,
#[Out]#          4.77315735e-11, 4.77312902e-11, 4.77310005e-11, 4.77307220e-11,
#[Out]#          4.77304227e-11] sr>)
h2cssubcube_common_beam = h2cssubcube.beams.common_beam()
h2cssubcube_b = h2cssubcube.convolve_to(h2cssubcube_common_beam)
targetsubcube_common_beam = targetsubcube.beams.common_beam()
targetsubcube_b = targetsubcube.convolve_to(targetsubcube_common_beam)
h2cssubcube_b.beam, targetsubcube_b.beam
#[Out]# (Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg,
#[Out]#  Beam: BMAJ=1.510909080505371 arcsec BMIN=1.1864451169967651 arcsec BPA=-81.83074951171875 deg)
# This might not be the neatest way to do this
h2cssubcube_b = h2cssubcube_b.with_spectral_unit(u.km/u.s)
targetsubcube_b = targetsubcube_b.with_spectral_unit(u.km/u.s)
h2cssubcube_b, targetsubcube_b
#[Out]# (DaskSpectralCube with shape=(38, 512, 512) and unit=Jy / beam and chunk size (38, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     38  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(41, 512, 512) and unit=Jy / beam and chunk size (41, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     41  type_s: VRAD      unit_s: km / s  range:       -9.650 km / s:      30.136 km / s)
import numpy as np
velocity_res_1 = np.diff(h2cssubcube_b.spectral_axis)[0]
velocity_res_2 = np.diff(targetsubcube_b.spectral_axis)[0]
np.abs(velocity_res_1), np.abs(velocity_res_2)
#[Out]# (<Quantity 1.08196348 km / s>, <Quantity 0.99465058 km / s>)
fwhm_gaussian = (velocity_res_1**2 - velocity_res_2**2)**0.5
fwhm_gaussian
#[Out]# <Quantity 0.42581121 km / s>
from astropy.convolution import Gaussian1DKernel
fwhm_to_sigma = np.sqrt(8*np.log(2))
# We want the kernel in pixel units, so we force to km/s and take the value
spectral_smoothing_kernel = Gaussian1DKernel(stddev=fwhm_gaussian.to(u.km/u.s).value / fwhm_to_sigma)
vel_lo = np.max([h2cssubcube_b.spectral_axis.min().value, 
                 targetsubcube_b.spectral_axis.min().value])*u.km/u.s
vel_hi = np.min([h2cssubcube_b.spectral_axis.max().value, 
                 targetsubcube_b.spectral_axis.max().value])*u.km/u.s

h2cssubcube_bc = h2cssubcube_b.spectral_slab(vel_lo, vel_hi)
targetsubcube_bc = targetsubcube_b.spectral_slab(vel_lo, vel_hi)
h2cssubcube_bc, targetsubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(40, 512, 512) and unit=Jy / beam and chunk size (40, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     40  type_s: VRAD      unit_s: km / s  range:       -9.650 km / s:      29.141 km / s)
targetsubcube_bc_spec = targetsubcube_bc.spectral_smooth(spectral_smoothing_kernel)
targetsubcube_bc_spec_resample = targetsubcube_bc_spec.spectral_interpolate(h2cssubcube_bc.spectral_axis)
targetsubcube_bc_spec_resample
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
h2cssubcube_bc.beam, targetsubcube_bc_spec_resample.beam
#[Out]# (Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg,
#[Out]#  Beam: BMAJ=1.510909080505371 arcsec BMIN=1.1864451169967651 arcsec BPA=-81.83074951171875 deg)
import radio_beam
common_beam = radio_beam.commonbeam.common_2beams(radio_beam.Beams(beams=[h2cssubcube_bc.beam, 
                                                                          targetsubcube_bc_spec_resample.beam]))
common_beam
#[Out]# Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg
targetsubcube_bc_spec_resample
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
# h2cssubcube_bc = h2cssubcube_bc.to(u.K) # Is this allowed???
targetsubcube_bc_spec_resample_spat = targetsubcube_bc_spec_resample.to(u.K).convolve_to(common_beam)
targetsubcube_bc_spec_resample_spat
# This takes a long time, and has a new warning:
# WARNING: nan_treatment='interpolate', however, NaN values detected post convolution.
# A contiguous region of NaN values, larger than the kernel size, are present in the input array. 
# Increase the kernel size to avoid this. [astropy.convolution.convolve]
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
import astropy
print(astropy.__version__)
import reproject
print(reproject.__version__)
# Need development version of astropy
# 4.3.dev1788+ga3263b6 works for me
