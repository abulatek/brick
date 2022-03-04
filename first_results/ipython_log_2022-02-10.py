########################################################
# Started Logging At: 2022-02-10 14:30:22
########################################################
########################################################
# # Started Logging At: 2022-02-10 14:30:23
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
h2cscube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 133582251086.390 Hz:135456817280.560 Hz
# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
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
import astropy
print(astropy.__version__)
import reproject
print(reproject.__version__)
# Need development version of astropy
# 4.3.dev1788+ga3263b6 works for me
targetsubcube_bc_reproj = targetsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
targetsubcube_bc_reproj
# This is currently taking a longgggg time to run, which is a sign that something went wrong?
end = time.time()
print(f"Time elapsed: {end - start} seconds")
import astropy
print(astropy.__version__)
import reproject
print(reproject.__version__)
import spectral_cube
print(spectral_cube.__version__)
import dask
print(dask.__version__)
# Need development version of astropy
# 4.3.dev1788+ga3263b6 works for me

# AG is:
#5.1.dev541+ged8cab8
#0.9.dev20+g41dbdf3
targetsubcube_bc_reproj, h2cssubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s)
h2cssubcube_bc[:,256,256].to(u.K).with_spectral_unit(u.km/u.s).quicklook()
targetsubcube_bc_reproj[:,256,256].quicklook()
med1 = h2cssubcube_bc.median(axis=0)  
h2cssubcube_f = h2cssubcube_bc - med1
med2 = targetsubcube_bc_reproj.median(axis=0)  
targetsubcube_f_reproj = targetsubcube_bc_reproj - med2
import matplotlib.pyplot as plt

# For rectangle:
import matplotlib.patches as mpatches
fig = plt.figure(dpi = display_dpi)
ax = fig.add_subplot(111)

# Make signal mask out of template molecule cube, in noise-free area (should be flat)
h2cs_sclip = h2cssubcube_f.sigma_clip_spectrally(3)
mad_std_spectrum_sclip = h2cs_sclip.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
         drawstyle = 'steps-mid', c = 'k')
plt.ylim(0.00175, 0.00245)
plt.xlabel(f'Velocity [{mad_std_spectrum_sclip.spectral_axis.unit.to_string("latex_inline")}]')
plt.ylabel(f'Standard deviation [{mad_std_spectrum_sclip.unit.to_string("latex_inline")}]')

flat_noise = mpatches.Rectangle((-7, 0.0017), 14, 0.001, alpha = 0.2, facecolor = "tab:blue")
plt.gca().add_patch(flat_noise)

plt.savefig("figures/flat_noise.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
h2cs_sclip_cut = h2cs_sclip.spectral_slab(-7*u.km/u.s, 7*u.km/u.s)
mad_std_spectrum_sclip_cut = h2cs_sclip_cut.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip_cut.spectral_axis.value, mad_std_spectrum_sclip_cut.value, 
         drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r'Standard deviation (K)')
mad_std_map_sclip = h2cs_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
mad_std_map_sclip.write('temperature_map/methyl_cyanide/template_noise.fits', overwrite = True)
mad_std_map_sclip.quicklook()
plain_mask = h2cssubcube_f >= 3 * mad_std_map_sclip # Get plain 3sigma mask
plain_masked_slab = h2cssubcube_f.with_mask(plain_mask) # Mask the template molecule cube
import scipy.ndimage as nd
# Make a low and high mask
low_snr_mask = (plain_masked_slab > 3 * mad_std_map_sclip).include()
high_snr_mask = (plain_masked_slab > 10 * mad_std_map_sclip).include()
low_snr_mask = low_snr_mask.compute() # Don't need this for this tutorial
high_snr_mask = high_snr_mask.compute() # Don't need this for this tutorial
# Find connected structures
structure = np.ones((3, 3, 3), dtype=bool)
low_snr_mask_labels, num_labels = nd.label(low_snr_mask, structure=structure)
print(f"Initial number of regions found: {num_labels}")
# From the labels, count the number of pixels within each label.
num_pixels_in_high_snr_mask = nd.sum(high_snr_mask,
                                     labels=low_snr_mask_labels,
                                     index=range(1, num_labels + 1)) # +1 offset for mask labels
# Repeat for the high signal mask.
num_pixels_in_low_snr_mask = nd.sum(low_snr_mask,
                                    labels=low_snr_mask_labels,
                                    index=range(1, num_labels + 1)) # +1 offset for mask labels
# To preserve the low_snr_mask, we will create a new signal mask where we will remove 
# regions that do not pass the criteria.
signal_mask = low_snr_mask
low_min_pixels = 40
high_min_pixels = 10
for num, (high_pix_num, low_pix_num) in enumerate(zip(num_pixels_in_high_snr_mask, 
                                                      num_pixels_in_low_snr_mask)):
    if high_pix_num >= high_min_pixels and low_pix_num >= low_min_pixels:
        # This region passes the criteria. Keep it in the mask.
        continue
    # Remove regions that do not pass the criteria.
    # NOTE: enumerate will start with 0, but the mask labels start at 1
    # We apply a +1 offset to `num` to account for this.
    signal_mask[low_snr_mask_labels == num + 1] = False
signal_mask_labels, num_labels = nd.label(signal_mask,
                                          structure=structure)
print(f"Final number of regions found: {num_labels}")
signal_mask = nd.binary_dilation(signal_mask, structure=structure, iterations=1)
masked_cube = targetsubcube_f_reproj.with_mask(signal_mask)
# We can write the masked cube to a file and look at it in DS9, but this is optional:
# masked_cube.write('masked_target_cube.fits', overwrite=True)
# masked_cube = SpectralCube.read('masked_target_cube.fits')
# Imports and matplotlib settings
from astropy.visualization import simple_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['figure.facecolor'] = 'w'

# Collapse the mask in three dimensions
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

# Set up maximum value for colorbar
vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)
cm = pl.matplotlib.cm.inferno.copy()
cm.set_under('w') # Make sure the "zero" color is white

fig = plt.figure(figsize = (8.5, 10)) # constrained_layout=True

# Get WCS coordinates from an (arbitrary?) masked cube
masked_moment0 = masked_cube.moment0()
wcs = masked_moment0.wcs

# Collapse along axis 0
ax1 = plt.subplot(323, projection = wcs, aspect = 1)
ax1.set_aspect(1)
im1 = ax1.imshow(mask, origin = 'lower', cmap=cm, vmin=0.1)
ax1.tick_params(direction = 'in')
ax1.coords[0].set_major_formatter('hh:mm:ss')
ax1.coords[1].set_major_formatter('d.dd')
ax1.set_xlabel('Right ascension'), ax1.set_ylabel('Declination')
cbar1 = plt.colorbar(mappable = im1, ax = ax1)
cbar1.set_label('# of valid pixels along axis')

# This may not be right:
# overlay = ax1.get_coords_overlay('galactic')
# overlay.grid(color='white', ls='dotted')
# overlay[0].set_axislabel('Galactic longitude')
# overlay[1].set_axislabel('Galactic latitude')

# Collapse along axis 1
ax2 = plt.subplot(321, sharex = ax1) # , adjustable='box'
im2 = ax2.imshow(collapse1, origin = 'lower', norm=norm, cmap=cm)
ax2.xaxis.set_tick_params(labelbottom = False) # Remove redundant tick labels
ax2.tick_params(direction = 'in')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
ax2.set_ylabel('Velocity (km/s)')
cbar2 = plt.colorbar(mappable=im2, ax=ax2)
cbar2.set_label('# of valid pixels along axis')

# Collapse along axis 2
ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin = 'lower', norm=norm, cmap=cm)
ax3.yaxis.set_tick_params(labelleft = False) # Remove redundant tick labels
ax3.tick_params(direction = 'in')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])
ax3.set_xlabel('Velocity (km/s)')
cbar3 = plt.colorbar(mappable=im3, ax=ax3)
cbar3.set_label('# of valid pixels along axis')

plt.savefig("figures/mask.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
masked_moment0 = masked_cube.moment0()

ax = plt.subplot(projection=masked_moment0.wcs)
im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label('Integrated Intensity (K km/s)')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
masked_moment1 = masked_cube.moment1()

ax = plt.subplot(projection=masked_moment1.wcs)
im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
cbar = plt.colorbar(im)
cbar.set_label('Centroid (km/s)')

ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
from pylab import imshow
v_thresh = 1000
masked_moment1 = masked_cube.moment1()
masked_moment1_outliers = (masked_moment1 > v_thresh*u.km/u.s)|(masked_moment1 < -v_thresh*u.km/u.s)
imshow(masked_moment1_outliers, origin='lower') 
# Clumps of outliers might mean they're real, just outside of vel range
#[Out]# <matplotlib.image.AxesImage at 0x2afce2357d30>
max_vel_coord = np.unravel_index(np.nanargmin(masked_moment1), masked_moment1.shape)
spectrum = masked_cube[:, max_vel_coord[0], max_vel_coord[1]]
print(masked_moment1[max_vel_coord[0], max_vel_coord[1]])
plt.plot(spectrum.spectral_axis, spectrum.value, drawstyle='steps-mid')
plt.xlabel('Velocity')
plt.ylabel('Intensity')
#[Out]# Text(0, 0.5, 'Intensity')
mom0 = masked_cube.moment0()
mom0_mask = mom0 > 0.3*u.K*u.km/u.s # Mask pixels with mom0 less than threshold
print(f"Found {mom0_mask.sum()} good pixels")
masked_cube_no_outliers = masked_cube.with_mask(mom0_mask)
masked_cube_no_outliers
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
from astropy.io import fits

data = masked_cube_no_outliers # 'DaskSpectralCube' object has no attribute 'value'
header = fits.PrimaryHDU()
header['BUNIT'] = masked_cube_no_outliers.unit.to_string('fits')
fits.PrimaryHDU(data=data, header=header).writeto('temperature_map/methyl_cyanide/ch3cn_0_masked.fits') # Export cube
########################################################
# Started Logging At: 2022-02-10 14:33:32
########################################################
########################################################
# # Started Logging At: 2022-02-10 14:33:33
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
h2cscube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 133582251086.390 Hz:135456817280.560 Hz
# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
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
import spectral_cube
print(spectral_cube.__version__)
import dask
print(dask.__version__)
# Need development version of astropy
# 4.3.dev1788+ga3263b6 works for me

# AG is:
#5.1.dev541+ged8cab8
#0.9.dev20+g41dbdf3
#0.6.1.dev38+g52cfdb9
#2021.10.0
start2 = time.time()
targetsubcube_bc_reproj = targetsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
targetsubcube_bc_reproj
# This is currently taking a longgggg time to run, which is a sign that something went wrong?
end = time.time()
print(f"Time elapsed: {end - start2} seconds (total since start = {end-start} seconds)")
targetsubcube_bc_reproj, h2cssubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s)
h2cssubcube_bc[:,256,256].to(u.K).with_spectral_unit(u.km/u.s).quicklook()
targetsubcube_bc_reproj[:,256,256].quicklook()
med1 = h2cssubcube_bc.median(axis=0)  
h2cssubcube_f = h2cssubcube_bc - med1
med2 = targetsubcube_bc_reproj.median(axis=0)  
targetsubcube_f_reproj = targetsubcube_bc_reproj - med2
import matplotlib.pyplot as plt

# For rectangle:
import matplotlib.patches as mpatches
fig = plt.figure(dpi = display_dpi)
ax = fig.add_subplot(111)

# Make signal mask out of template molecule cube, in noise-free area (should be flat)
h2cs_sclip = h2cssubcube_f.sigma_clip_spectrally(3)
mad_std_spectrum_sclip = h2cs_sclip.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
         drawstyle = 'steps-mid', c = 'k')
plt.ylim(0.00175, 0.00245)
plt.xlabel(f'Velocity [{mad_std_spectrum_sclip.spectral_axis.unit.to_string("latex_inline")}]')
plt.ylabel(f'Standard deviation [{mad_std_spectrum_sclip.unit.to_string("latex_inline")}]')

flat_noise = mpatches.Rectangle((-7, 0.0017), 14, 0.001, alpha = 0.2, facecolor = "tab:blue")
plt.gca().add_patch(flat_noise)

plt.savefig("figures/flat_noise.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
h2cs_sclip_cut = h2cs_sclip.spectral_slab(-7*u.km/u.s, 7*u.km/u.s)
mad_std_spectrum_sclip_cut = h2cs_sclip_cut.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip_cut.spectral_axis.value, mad_std_spectrum_sclip_cut.value, 
         drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r'Standard deviation (K)')
mad_std_map_sclip = h2cs_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
mad_std_map_sclip.write('temperature_map/methyl_cyanide/template_noise.fits', overwrite = True)
mad_std_map_sclip.quicklook()
plain_mask = h2cssubcube_f >= 3 * mad_std_map_sclip # Get plain 3sigma mask
plain_masked_slab = h2cssubcube_f.with_mask(plain_mask) # Mask the template molecule cube
import scipy.ndimage as nd
# Make a low and high mask
low_snr_mask = (plain_masked_slab > 3 * mad_std_map_sclip).include()
high_snr_mask = (plain_masked_slab > 10 * mad_std_map_sclip).include()
low_snr_mask = low_snr_mask.compute() # Don't need this for this tutorial
high_snr_mask = high_snr_mask.compute() # Don't need this for this tutorial
# Find connected structures
structure = np.ones((3, 3, 3), dtype=bool)
low_snr_mask_labels, num_labels = nd.label(low_snr_mask, structure=structure)
print(f"Initial number of regions found: {num_labels}")
# From the labels, count the number of pixels within each label.
num_pixels_in_high_snr_mask = nd.sum(high_snr_mask,
                                     labels=low_snr_mask_labels,
                                     index=range(1, num_labels + 1)) # +1 offset for mask labels
# Repeat for the high signal mask.
num_pixels_in_low_snr_mask = nd.sum(low_snr_mask,
                                    labels=low_snr_mask_labels,
                                    index=range(1, num_labels + 1)) # +1 offset for mask labels
# To preserve the low_snr_mask, we will create a new signal mask where we will remove 
# regions that do not pass the criteria.
signal_mask = low_snr_mask
low_min_pixels = 40
high_min_pixels = 10
for num, (high_pix_num, low_pix_num) in enumerate(zip(num_pixels_in_high_snr_mask, 
                                                      num_pixels_in_low_snr_mask)):
    if high_pix_num >= high_min_pixels and low_pix_num >= low_min_pixels:
        # This region passes the criteria. Keep it in the mask.
        continue
    # Remove regions that do not pass the criteria.
    # NOTE: enumerate will start with 0, but the mask labels start at 1
    # We apply a +1 offset to `num` to account for this.
    signal_mask[low_snr_mask_labels == num + 1] = False
signal_mask_labels, num_labels = nd.label(signal_mask,
                                          structure=structure)
print(f"Final number of regions found: {num_labels}")
signal_mask = nd.binary_dilation(signal_mask, structure=structure, iterations=1)
masked_cube = targetsubcube_f_reproj.with_mask(signal_mask)
# We can write the masked cube to a file and look at it in DS9, but this is optional:
# masked_cube.write('masked_target_cube.fits', overwrite=True)
# masked_cube = SpectralCube.read('masked_target_cube.fits')
# Imports and matplotlib settings
from astropy.visualization import simple_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['figure.facecolor'] = 'w'

# Collapse the mask in three dimensions
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

# Set up maximum value for colorbar
vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)
cm = pl.matplotlib.cm.inferno.copy()
cm.set_under('w') # Make sure the "zero" color is white

fig = plt.figure(figsize = (8.5, 10)) # constrained_layout=True

# Get WCS coordinates from an (arbitrary?) masked cube
masked_moment0 = masked_cube.moment0()
wcs = masked_moment0.wcs

# Collapse along axis 0
ax1 = plt.subplot(323, projection = wcs, aspect = 1)
ax1.set_aspect(1)
im1 = ax1.imshow(mask, origin = 'lower', cmap=cm, vmin=0.1)
ax1.tick_params(direction = 'in')
ax1.coords[0].set_major_formatter('hh:mm:ss')
ax1.coords[1].set_major_formatter('d.dd')
ax1.set_xlabel('Right ascension'), ax1.set_ylabel('Declination')
cbar1 = plt.colorbar(mappable = im1, ax = ax1)
cbar1.set_label('# of valid pixels along axis')

# This may not be right:
# overlay = ax1.get_coords_overlay('galactic')
# overlay.grid(color='white', ls='dotted')
# overlay[0].set_axislabel('Galactic longitude')
# overlay[1].set_axislabel('Galactic latitude')

# Collapse along axis 1
ax2 = plt.subplot(321, sharex = ax1) # , adjustable='box'
im2 = ax2.imshow(collapse1, origin = 'lower', norm=norm, cmap=cm)
ax2.xaxis.set_tick_params(labelbottom = False) # Remove redundant tick labels
ax2.tick_params(direction = 'in')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
ax2.set_ylabel('Velocity (km/s)')
cbar2 = plt.colorbar(mappable=im2, ax=ax2)
cbar2.set_label('# of valid pixels along axis')

# Collapse along axis 2
ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin = 'lower', norm=norm, cmap=cm)
ax3.yaxis.set_tick_params(labelleft = False) # Remove redundant tick labels
ax3.tick_params(direction = 'in')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])
ax3.set_xlabel('Velocity (km/s)')
cbar3 = plt.colorbar(mappable=im3, ax=ax3)
cbar3.set_label('# of valid pixels along axis')

plt.savefig("figures/mask.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
masked_moment0 = masked_cube.moment0()

ax = plt.subplot(projection=masked_moment0.wcs)
im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label('Integrated Intensity (K km/s)')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
masked_moment1 = masked_cube.moment1()

ax = plt.subplot(projection=masked_moment1.wcs)
im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
cbar = plt.colorbar(im)
cbar.set_label('Centroid (km/s)')

ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
from pylab import imshow
v_thresh = 1000
masked_moment1 = masked_cube.moment1()
masked_moment1_outliers = (masked_moment1 > v_thresh*u.km/u.s)|(masked_moment1 < -v_thresh*u.km/u.s)
imshow(masked_moment1_outliers, origin='lower') 
# Clumps of outliers might mean they're real, just outside of vel range
#[Out]# <matplotlib.image.AxesImage at 0x2add0eaf0070>
max_vel_coord = np.unravel_index(np.nanargmin(masked_moment1), masked_moment1.shape)
spectrum = masked_cube[:, max_vel_coord[0], max_vel_coord[1]]
print(masked_moment1[max_vel_coord[0], max_vel_coord[1]])
plt.plot(spectrum.spectral_axis, spectrum.value, drawstyle='steps-mid')
plt.xlabel('Velocity')
plt.ylabel('Intensity')
#[Out]# Text(0, 0.5, 'Intensity')
mom0 = masked_cube.moment0()
mom0_mask = mom0 > 0.3*u.K*u.km/u.s # Mask pixels with mom0 less than threshold
print(f"Found {mom0_mask.sum()} good pixels")
masked_cube_no_outliers = masked_cube.with_mask(mom0_mask)
masked_cube_no_outliers
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
from astropy.io import fits

data = masked_cube_no_outliers # 'DaskSpectralCube' object has no attribute 'value'
header = fits.PrimaryHDU()
header['BUNIT'] = masked_cube_no_outliers.unit.to_string('fits')
fits.PrimaryHDU(data=data, header=header).writeto('temperature_map/methyl_cyanide/ch3cn_0_masked.fits') # Export cube
########################################################
# Started Logging At: 2022-02-10 14:35:22
########################################################
########################################################
# # Started Logging At: 2022-02-10 14:35:23
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
h2cscube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 133582251086.390 Hz:135456817280.560 Hz
# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
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
h2cssubcube = h2cssubcube.to(u.K)
targetsubcube = targetsubcube.to(u.K)
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
#[Out]# (DaskSpectralCube with shape=(38, 512, 512) and unit=K and chunk size (38, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     38  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(41, 512, 512) and unit=K and chunk size (41, 512, 512):
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
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(40, 512, 512) and unit=K and chunk size (40, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     40  type_s: VRAD      unit_s: km / s  range:       -9.650 km / s:      29.141 km / s)
targetsubcube_bc_spec = targetsubcube_bc.spectral_smooth(spectral_smoothing_kernel)
targetsubcube_bc_spec_resample = targetsubcube_bc_spec.spectral_interpolate(h2cssubcube_bc.spectral_axis)
targetsubcube_bc_spec_resample
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
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
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
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
import spectral_cube
print(spectral_cube.__version__)
import dask
print(dask.__version__)
# Need development version of astropy
# 4.3.dev1788+ga3263b6 works for me

# AG is:
#5.1.dev541+ged8cab8
#0.9.dev20+g41dbdf3
#0.6.1.dev38+g52cfdb9
#2021.10.0
start2 = time.time()
targetsubcube_bc_reproj = targetsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
targetsubcube_bc_reproj
# This is currently taking a longgggg time to run, which is a sign that something went wrong?
end = time.time()
print(f"Time elapsed: {end - start2} seconds (total since start = {end-start} seconds)")
targetsubcube_bc_reproj, h2cssubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s)
h2cssubcube_bc[:,256,256].to(u.K).with_spectral_unit(u.km/u.s).quicklook()
targetsubcube_bc_reproj[:,256,256].quicklook()
med1 = h2cssubcube_bc.median(axis=0)  
h2cssubcube_f = h2cssubcube_bc - med1
med2 = targetsubcube_bc_reproj.median(axis=0)  
targetsubcube_f_reproj = targetsubcube_bc_reproj - med2
import matplotlib.pyplot as plt

# For rectangle:
import matplotlib.patches as mpatches
fig = plt.figure(dpi = display_dpi)
ax = fig.add_subplot(111)

# Make signal mask out of template molecule cube, in noise-free area (should be flat)
h2cs_sclip = h2cssubcube_f.sigma_clip_spectrally(3)
mad_std_spectrum_sclip = h2cs_sclip.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
         drawstyle = 'steps-mid', c = 'k')
plt.ylim(0.00175, 0.00245)
plt.xlabel(f'Velocity [{mad_std_spectrum_sclip.spectral_axis.unit.to_string("latex_inline")}]')
plt.ylabel(f'Standard deviation [{mad_std_spectrum_sclip.unit.to_string("latex_inline")}]')

flat_noise = mpatches.Rectangle((-7, 0.0017), 14, 0.001, alpha = 0.2, facecolor = "tab:blue")
plt.gca().add_patch(flat_noise)

plt.savefig("figures/flat_noise.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
h2cs_sclip_cut = h2cs_sclip.spectral_slab(-7*u.km/u.s, 7*u.km/u.s)
mad_std_spectrum_sclip_cut = h2cs_sclip_cut.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip_cut.spectral_axis.value, mad_std_spectrum_sclip_cut.value, 
         drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r'Standard deviation (K)')
mad_std_map_sclip = h2cs_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
mad_std_map_sclip.write('temperature_map/methyl_cyanide/template_noise.fits', overwrite = True)
mad_std_map_sclip.quicklook()
plain_mask = h2cssubcube_f >= 3 * mad_std_map_sclip # Get plain 3sigma mask
plain_masked_slab = h2cssubcube_f.with_mask(plain_mask) # Mask the template molecule cube
import scipy.ndimage as nd
# Make a low and high mask
low_snr_mask = (plain_masked_slab > 3 * mad_std_map_sclip).include()
high_snr_mask = (plain_masked_slab > 10 * mad_std_map_sclip).include()
low_snr_mask = low_snr_mask.compute() # Don't need this for this tutorial
high_snr_mask = high_snr_mask.compute() # Don't need this for this tutorial
# Find connected structures
structure = np.ones((3, 3, 3), dtype=bool)
low_snr_mask_labels, num_labels = nd.label(low_snr_mask, structure=structure)
print(f"Initial number of regions found: {num_labels}")
# From the labels, count the number of pixels within each label.
num_pixels_in_high_snr_mask = nd.sum(high_snr_mask,
                                     labels=low_snr_mask_labels,
                                     index=range(1, num_labels + 1)) # +1 offset for mask labels
# Repeat for the high signal mask.
num_pixels_in_low_snr_mask = nd.sum(low_snr_mask,
                                    labels=low_snr_mask_labels,
                                    index=range(1, num_labels + 1)) # +1 offset for mask labels
# To preserve the low_snr_mask, we will create a new signal mask where we will remove 
# regions that do not pass the criteria.
signal_mask = low_snr_mask
low_min_pixels = 40
high_min_pixels = 10
for num, (high_pix_num, low_pix_num) in enumerate(zip(num_pixels_in_high_snr_mask, 
                                                      num_pixels_in_low_snr_mask)):
    if high_pix_num >= high_min_pixels and low_pix_num >= low_min_pixels:
        # This region passes the criteria. Keep it in the mask.
        continue
    # Remove regions that do not pass the criteria.
    # NOTE: enumerate will start with 0, but the mask labels start at 1
    # We apply a +1 offset to `num` to account for this.
    signal_mask[low_snr_mask_labels == num + 1] = False
signal_mask_labels, num_labels = nd.label(signal_mask,
                                          structure=structure)
print(f"Final number of regions found: {num_labels}")
signal_mask = nd.binary_dilation(signal_mask, structure=structure, iterations=1)
masked_cube = targetsubcube_f_reproj.with_mask(signal_mask)
# We can write the masked cube to a file and look at it in DS9, but this is optional:
# masked_cube.write('masked_target_cube.fits', overwrite=True)
# masked_cube = SpectralCube.read('masked_target_cube.fits')
# Imports and matplotlib settings
from astropy.visualization import simple_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['figure.facecolor'] = 'w'

# Collapse the mask in three dimensions
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

# Set up maximum value for colorbar
vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)
cm = pl.matplotlib.cm.inferno.copy()
cm.set_under('w') # Make sure the "zero" color is white

fig = plt.figure(figsize = (8.5, 10)) # constrained_layout=True

# Get WCS coordinates from an (arbitrary?) masked cube
masked_moment0 = masked_cube.moment0()
wcs = masked_moment0.wcs

# Collapse along axis 0
ax1 = plt.subplot(323, projection = wcs, aspect = 1)
ax1.set_aspect(1)
im1 = ax1.imshow(mask, origin = 'lower', cmap=cm, vmin=0.1)
ax1.tick_params(direction = 'in')
ax1.coords[0].set_major_formatter('hh:mm:ss')
ax1.coords[1].set_major_formatter('d.dd')
ax1.set_xlabel('Right ascension'), ax1.set_ylabel('Declination')
cbar1 = plt.colorbar(mappable = im1, ax = ax1)
cbar1.set_label('# of valid pixels along axis')

# This may not be right:
# overlay = ax1.get_coords_overlay('galactic')
# overlay.grid(color='white', ls='dotted')
# overlay[0].set_axislabel('Galactic longitude')
# overlay[1].set_axislabel('Galactic latitude')

# Collapse along axis 1
ax2 = plt.subplot(321, sharex = ax1) # , adjustable='box'
im2 = ax2.imshow(collapse1, origin = 'lower', norm=norm, cmap=cm)
ax2.xaxis.set_tick_params(labelbottom = False) # Remove redundant tick labels
ax2.tick_params(direction = 'in')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
ax2.set_ylabel('Velocity (km/s)')
cbar2 = plt.colorbar(mappable=im2, ax=ax2)
cbar2.set_label('# of valid pixels along axis')

# Collapse along axis 2
ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin = 'lower', norm=norm, cmap=cm)
ax3.yaxis.set_tick_params(labelleft = False) # Remove redundant tick labels
ax3.tick_params(direction = 'in')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])
ax3.set_xlabel('Velocity (km/s)')
cbar3 = plt.colorbar(mappable=im3, ax=ax3)
cbar3.set_label('# of valid pixels along axis')

plt.savefig("figures/mask.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
masked_moment0 = masked_cube.moment0()

ax = plt.subplot(projection=masked_moment0.wcs)
im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label('Integrated Intensity (K km/s)')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
masked_moment1 = masked_cube.moment1()

ax = plt.subplot(projection=masked_moment1.wcs)
im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
cbar = plt.colorbar(im)
cbar.set_label('Centroid (km/s)')

ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
from pylab import imshow
v_thresh = 1000
masked_moment1 = masked_cube.moment1()
masked_moment1_outliers = (masked_moment1 > v_thresh*u.km/u.s)|(masked_moment1 < -v_thresh*u.km/u.s)
imshow(masked_moment1_outliers, origin='lower') 
# Clumps of outliers might mean they're real, just outside of vel range
#[Out]# <matplotlib.image.AxesImage at 0x2ac30feb1b50>
max_vel_coord = np.unravel_index(np.nanargmin(masked_moment1), masked_moment1.shape)
spectrum = masked_cube[:, max_vel_coord[0], max_vel_coord[1]]
print(masked_moment1[max_vel_coord[0], max_vel_coord[1]])
plt.plot(spectrum.spectral_axis, spectrum.value, drawstyle='steps-mid')
plt.xlabel('Velocity')
plt.ylabel('Intensity')
#[Out]# Text(0, 0.5, 'Intensity')
mom0 = masked_cube.moment0()
mom0_mask = mom0 > 0.3*u.K*u.km/u.s # Mask pixels with mom0 less than threshold
print(f"Found {mom0_mask.sum()} good pixels")
masked_cube_no_outliers = masked_cube.with_mask(mom0_mask)
masked_cube_no_outliers
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
from astropy.io import fits

data = masked_cube_no_outliers # 'DaskSpectralCube' object has no attribute 'value'
header = fits.PrimaryHDU()
header['BUNIT'] = masked_cube_no_outliers.unit.to_string('fits')
fits.PrimaryHDU(data=data, header=header).writeto('temperature_map/methyl_cyanide/ch3cn_0_masked.fits') # Export cube
from dask.diagnostics import ProgressBar
ProgressBar().register()
start2 = time.time()
targetsubcube_bc_reproj = targetsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
targetsubcube_bc_reproj
# This is currently taking a longgggg time to run, which is a sign that something went wrong?
end = time.time()
print(f"Time elapsed: {end - start2} seconds (total since start = {end-start} seconds)")
