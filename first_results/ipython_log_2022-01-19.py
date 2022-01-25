########################################################
# Started Logging At: 2022-01-19 16:08:50
########################################################
########################################################
# # Started Logging At: 2022-01-19 16:08:50
########################################################
import pylab as pl
# pl.style.use('dark_background')

# plt.rc('text', usetex = True) # Use LaTeX font in plots
# plt.rc('font', family = 'serif')
# plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
# plt.rcParams.update({'font.size': 14})
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
from astroquery.splatalogue import Splatalogue
h2cscube.find_lines(chemical_name='H2CS', line_lists=['JPL'], 
                    show_upper_degeneracy=True).show_in_notebook()
########################################################
# Started Logging At: 2022-01-19 16:09:04
########################################################
########################################################
# # Started Logging At: 2022-01-19 16:09:05
########################################################
########################################################
# Started Logging At: 2022-01-19 16:09:06
########################################################
########################################################
# # Started Logging At: 2022-01-19 16:09:06
########################################################
import pylab as pl
# pl.style.use('dark_background')

# plt.rc('text', usetex = True) # Use LaTeX font in plots
# plt.rc('font', family = 'serif')
# plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
# plt.rcParams.update({'font.size': 14})
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
ch3cncube = ch3cncube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                         rest_value=ch3cn_freqs[0]*u.GHz).spectral_slab(-10*u.km/u.s, 
                                                                                        80*u.km/u.s)
print(ch3cncube)
ch3cnsubcube = ch3cncube.spectral_slab(-10*u.km/u.s, 30*u.km/u.s)
h2cssubcube.max(axis=0).quicklook()
# This first quicklook has to be re-run, as it doesn't run the first time (not sure why)
ch3cnsubcube.max(axis=0).quicklook()
h2cssubcube.beams, ch3cnsubcube.beams
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
ch3cnsubcube_common_beam = ch3cnsubcube.beams.common_beam()
ch3cnsubcube_b = ch3cnsubcube.convolve_to(ch3cnsubcube_common_beam)
h2cssubcube_b.beam, ch3cnsubcube_b.beam
#[Out]# (Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg,
#[Out]#  Beam: BMAJ=1.510909080505371 arcsec BMIN=1.1864451169967651 arcsec BPA=-81.83074951171875 deg)
# This might not be the neatest way to do this
h2cssubcube_b = h2cssubcube_b.with_spectral_unit(u.km/u.s)
ch3cnsubcube_b = ch3cnsubcube_b.with_spectral_unit(u.km/u.s)
h2cssubcube_b, ch3cnsubcube_b
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
velocity_res_2 = np.diff(ch3cnsubcube_b.spectral_axis)[0]
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
                 ch3cnsubcube_b.spectral_axis.min().value])*u.km/u.s
vel_hi = np.min([h2cssubcube_b.spectral_axis.max().value, 
                 ch3cnsubcube_b.spectral_axis.max().value])*u.km/u.s

h2cssubcube_bc = h2cssubcube_b.spectral_slab(vel_lo, vel_hi)
ch3cnsubcube_bc = ch3cnsubcube_b.spectral_slab(vel_lo, vel_hi)
h2cssubcube_bc, ch3cnsubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(40, 512, 512) and unit=Jy / beam and chunk size (40, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     40  type_s: VRAD      unit_s: km / s  range:       -9.650 km / s:      29.141 km / s)
ch3cnsubcube_bc_spec = ch3cnsubcube_bc.spectral_smooth(spectral_smoothing_kernel)
ch3cnsubcube_bc_spec_resample = ch3cnsubcube_bc_spec.spectral_interpolate(h2cssubcube_bc.spectral_axis)
ch3cnsubcube_bc_spec_resample
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
h2cssubcube_bc.beam, ch3cnsubcube_bc_spec_resample.beam
#[Out]# (Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg,
#[Out]#  Beam: BMAJ=1.510909080505371 arcsec BMIN=1.1864451169967651 arcsec BPA=-81.83074951171875 deg)
import radio_beam
common_beam = radio_beam.commonbeam.common_2beams(radio_beam.Beams(beams=[h2cssubcube_bc.beam, 
                                                                          ch3cnsubcube_bc_spec_resample.beam]))
common_beam
#[Out]# Beam: BMAJ=1.621724247932434 arcsec BMIN=1.2848083972930908 arcsec BPA=-84.60916900634766 deg
ch3cnsubcube_bc_spec_resample_spat = ch3cnsubcube_bc_spec_resample.to(u.K).convolve_to(common_beam)
ch3cnsubcube_bc_spec_resample_spat
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
ch3cnsubcube_bc_reproj = ch3cnsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
ch3cnsubcube_bc_reproj
# This is currently taking a longgggg time to run, which is a sign that something went wrong?
#[Out]# DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s
ch3cnsubcube_bc_reproj, h2cssubcube_bc
#[Out]# (DaskSpectralCube with shape=(37, 512, 512) and unit=K and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s,
#[Out]#  DaskSpectralCube with shape=(37, 512, 512) and unit=Jy / beam and chunk size (37, 512, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:     37  type_s: VRAD      unit_s: km / s  range:       -9.343 km / s:      29.607 km / s)
h2cssubcube_bc[:,256,256].to(u.K).with_spectral_unit(u.km/u.s).quicklook()
ch3cnsubcube_bc_reproj[:,256,256].quicklook()
med1 = h2cssubcube_bc.median(axis=0)  
h2cssubcube_f = h2cssubcube_bc - med1
med2 = ch3cnsubcube_bc_reproj.median(axis=0)  
ch3cnsubcube_f_reproj = ch3cnsubcube_bc_reproj - med2
import matplotlib.pyplot as plt

# For rectangle:
import matplotlib.patches as mpatches
fig = plt.figure(dpi = display_dpi)
ax = fig.add_subplot(111)

# Make signal mask out of template molecule cube, in noise-free area (should be flat)
h2cs_sclip = h2cssubcube_f.sigma_clip_spectrally(3)
mad_std_spectrum_sclip = h2cs_sclip.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
         drawstyle='steps-mid', c='k')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Standard deviation (K)')
plt.tight_layout()

flat_noise = mpatches.Rectangle((-7, 0.0017), 14, 0.001, alpha = 0.2, facecolor = "tab:blue")
plt.gca().add_patch(flat_noise)

plt.savefig("figures/flat_noise.pdf", dpi = 300, facecolor='w', edgecolor='w')
plt.show()
h2cs_sclip_cut = h2cs_sclip.spectral_slab(-7*u.km/u.s, 7*u.km/u.s)
mad_std_spectrum_sclip_cut = h2cs_sclip_cut.mad_std(axis=(1, 2))
plt.plot(mad_std_spectrum_sclip_cut.spectral_axis.value, mad_std_spectrum_sclip_cut.value, 
         drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r'Standard deviation (K)')
mad_std_map_sclip = h2cs_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
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
masked_cube = ch3cnsubcube_f_reproj.with_mask(signal_mask)
# We can write the masked cube to a file and look at it in DS9, but this is optional:
# masked_cube.write('masked_target_cube.fits', overwrite=True)
# masked_cube = SpectralCube.read('masked_target_cube.fits')
mask = signal_mask.astype('int')[25]
collapse1 = signal_mask.mean(axis = 1)
collapse2 = signal_mask.mean(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box', aspect=100)
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
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
#[Out]# <matplotlib.image.AxesImage at 0x2bab36accd90>
max_vel_coord = np.unravel_index(np.nanargmin(masked_moment1), masked_moment1.shape)
spectrum = masked_cube[:, max_vel_coord[0], max_vel_coord[1]]
print(masked_moment1[max_vel_coord[0], max_vel_coord[1]])
plt.plot(spectrum.spectral_axis, spectrum.value, drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Intensity (K)')
#[Out]# Text(0, 0.5, 'Intensity (K)')
mom0 = masked_cube.moment0()
mom0_mask = mom0 > 0.3*u.K*u.km/u.s # Mask pixels with mom0 less than threshold
print(f"Found {mom0_mask.sum()} good pixels")
masked_cube_no_outliers = masked_cube.with_mask(mom0_mask)
masked_moment0 = masked_cube_no_outliers.moment0()

ax = plt.subplot(projection=masked_moment0.wcs)
im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label('Integrated Intensity (K km/s)')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
masked_moment1 = masked_cube_no_outliers.moment1()

ax = plt.subplot(projection=masked_moment1.wcs)
im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
cbar = plt.colorbar(im)
cbar.set_label('Centroid (km/s)')

ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
masked_cube_no_outliers.shape # (75, 250, 250)
collapse1 = masked_cube_no_outliers.sum(axis = 1) # Collapse along one spatial axis
collapse2 = masked_cube_no_outliers.sum(axis = 2) # Collapse along the other spatial axis
imshow(collapse1.value, origin='lower')
#[Out]# <matplotlib.image.AxesImage at 0x2bab1c8969a0>
imshow(collapse2.value, origin='lower')
#[Out]# <matplotlib.image.AxesImage at 0x2bab33fc0b80>
# Some leftover code from going over this with Adam
# Log file is located here: 
# /orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/ipython_log_2021-06-21.py

# subcube_regrid.with_mask(mask, tolerance=1000)
# subcube_regrid.with_mask(mask.include())
# mask.include()
# mask = (h2cscube > 5*u.mJy/u.beam).include().compute()
# subcube_regrid.with_mask(mask)
# #subcube_regrid.with_mask(mask).write('
# get_ipython().run_line_magic('history', '')
# subcube_regrid.with_mask(mask).write('CH3CN_8(0)-7(0)_masked.fits')
# masked_subcube_regrid = subcube_regrid.with_mask(mask)
# masked_subcube_regrid.write('CH3CN_8(0)-7(0)_masked.fits', overwrite=True)
# masked_subcube_regrid
# masked_subcube_regrid.rechunk()
# masked_subcube_regrid.rechunk().write('CH3CN_8(0)-7(0)_masked.fits', overwrite=True)

# ch3cnK0 = SpectralCube.read('./thioformaldehyde/H2CS_303-202.fits')
# ch3cnK0
# ch3cnK0_rg = ch3cnK0.convolve_to(h2cscubes.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# ch3cnK0_rg.with_mask(mask)
# ch3cnK0 = SpectralCube.read('CH3CN_8(0)-7(0).fits')
# ch3cnK0_rg = ch3cnK0.convolve_to(h2cscubes.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# h2cscube
# ch3cnK0_rg = ch3cnK0.convolve_to(h2cscube.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# ch3cnK0_rg
# ch3cnK0_rg.with_mask(mask).write('CH3CN_8(0)-7(0)_masked.fits', overwrite=True)
# #cube = SpectralCube.read('../imaging_results/source_ab_137_spw69_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image')
# ch3cnK0_dask_cube = SpectralCube.read('../imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image').with_spect
# ral_unit(u.km/u.s, velocity_convention='radio', rest_value=147.1745883*u.GHz)
# ch3cnK0_dask_cube = SpectralCube.read('../imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image').with_spect
# ral_unit(u.km/u.s, velocity_convention='radio', rest_value=147.1745883*u.GHz)
# ch3cnK0_dc_rg = ch3cnK0_dask_cube.convolve_to(h2cscubes.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# ch3cnK0_dask_cube = SpectralCube.read('../imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image').with_spect
# ral_unit(u.km/u.s, velocity_convention='radio', rest_value=147.1745883*u.GHz)
# ch3cnK0_dc_rg = ch3cnK0_dask_cube.convolve_to(h2cscube.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# # the above is very fast b/c it was lazy
# ch3cnK0_dask_cube = SpectralCube.read('../imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=147.1745883*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# ch3cnK0_dc_rg = ch3cnK0_dask_cube.convolve_to(h2cscube.beams.common_beam()).spectral_interpolate(h2cscube.spectral_axis)
# # the above is very fast b/c it was lazy
# ch3cnK0_dc_rg
# ch3cnK0_dc_rg._data
# mch3cnK0_dc_rg = ch3cnK0_dc_rg.with_mask(mask)
# with ProgressBar():
#     mch3cnK0_dc_rg.write('TEST.fits')
# template_noise = mad_std_map_sclip
# template_noise.write('methyl_cyanide/template_noise.fits')
# ch3cn_7_masked = masked_cube_no_outliers
# ch3cn_7_masked.write('methyl_cyanide/ch3cn_7_masked.fits')
# This is an attempt to re-run the notebook for all rungs except index 0. It is untested.
# for ch3cn_freq in ch3cn_freqs[1:]:
    
#     # Convert cube spectral axes from frequency to velocity
#     ch3cncube = ch3cncube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
#                                              rest_value=ch3cn_freq*u.GHz).spectral_slab(-10*u.km/u.s, 
#                                                                                         80*u.km/u.s)
#     print(ch3cncube)
#     ch3cnsubcube = ch3cncube.spectral_slab(-10*u.km/u.s, 30*u.km/u.s)
    
#     # Do some quicklooks of the peak intensity to see what we expect to see
#     ch3cnsubcube.max(axis=0).quicklook()
    
#     # Spectral interpolation
#     ch3cnsubcube_common_beam = ch3cnsubcube.beams.common_beam()
#     ch3cnsubcube_b = ch3cnsubcube.convolve_to(ch3cnsubcube_common_beam)
#     ch3cnsubcube_b = ch3cnsubcube_b.with_spectral_unit(u.km/u.s)
#     velocity_res_2 = np.diff(ch3cnsubcube_b.spectral_axis)[0]
#     fwhm_gaussian = (velocity_res_1**2 - velocity_res_2**2)**0.5
#     spectral_smoothing_kernel = Gaussian1DKernel(stddev=fwhm_gaussian.to(u.km/u.s).value / fwhm_to_sigma)
#     vel_lo = np.max([h2cssubcube_b.spectral_axis.min().value, 
#                      ch3cnsubcube_b.spectral_axis.min().value])*u.km/u.s
#     vel_hi = np.min([h2cssubcube_b.spectral_axis.max().value, 
#                      ch3cnsubcube_b.spectral_axis.max().value])*u.km/u.s
#     h2cssubcube_bc = h2cssubcube_b.spectral_slab(vel_lo, vel_hi)
#     ch3cnsubcube_bc = ch3cnsubcube_b.spectral_slab(vel_lo, vel_hi)
#     ch3cnsubcube_bc_spec = ch3cnsubcube_bc.spectral_smooth(spectral_smoothing_kernel)
#     ch3cnsubcube_bc_spec_resample = ch3cnsubcube_bc_spec.spectral_interpolate(h2cssubcube_bc.spectral_axis)

#     # Spatial smoothing
#     ## Here we assume that the H2CS beam is larger, which might not always be right?
#     common_beam = radio_beam.commonbeam.common_2beams(radio_beam.Beams(beams=[h2cssubcube_bc.beam, 
#                                                                           ch3cnsubcube_bc_spec_resample.beam]))
#     ch3cnsubcube_bc_spec_resample_spat = ch3cnsubcube_bc_spec_resample.to(u.K).convolve_to(common_beam)
    
#     # Reprojection
#     ch3cnsubcube_bc_reproj = ch3cnsubcube_bc_spec_resample_spat.reproject(h2cssubcube_bc.header)
    
#     # Continuum subtraction
#     med1 = h2cssubcube_bc.median(axis=0)  
#     h2cssubcube_f = h2cssubcube_bc - med1
#     med2 = ch3cnsubcube_bc_reproj.median(axis=0)  
#     ch3cnsubcube_f_reproj = ch3cnsubcube_bc_reproj - med2
    
#     # More to do here...
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box', aspect=100)
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
pl.rcParams['figure.facecolor'] = 'w'
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box', aspect=100)
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box', aspect=100)
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(100)

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(10)

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[0]/collapse1.shape[1])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse1.shape[1]/collapse1.shape[0])



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])



# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

pl.colorbar()


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

pl.colorbar(mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab36e78130>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)


from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask, collapse1, collapse2])

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower') # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower') # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower') # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', vmax=vmax) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', vmax=vmax) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', vmax=vmax) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
from astropy.visualization import simple_norm
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax)

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)


# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax)

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab378f3be0>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax)

cm = pl.colors.viridis()

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm, cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax)

cm = pl.colormaps.viridis()

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm, cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax)

cm = pl.matplotlib.colors.viridis()

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm, cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
pl.colormaps
#[Out]# <matplotlib.cm.ColormapRegistry at 0x2baa93385100>
pl.colormaps()
#[Out]# ['magma',
#[Out]#  'inferno',
#[Out]#  'plasma',
#[Out]#  'viridis',
#[Out]#  'cividis',
#[Out]#  'twilight',
#[Out]#  'twilight_shifted',
#[Out]#  'turbo',
#[Out]#  'Blues',
#[Out]#  'BrBG',
#[Out]#  'BuGn',
#[Out]#  'BuPu',
#[Out]#  'CMRmap',
#[Out]#  'GnBu',
#[Out]#  'Greens',
#[Out]#  'Greys',
#[Out]#  'OrRd',
#[Out]#  'Oranges',
#[Out]#  'PRGn',
#[Out]#  'PiYG',
#[Out]#  'PuBu',
#[Out]#  'PuBuGn',
#[Out]#  'PuOr',
#[Out]#  'PuRd',
#[Out]#  'Purples',
#[Out]#  'RdBu',
#[Out]#  'RdGy',
#[Out]#  'RdPu',
#[Out]#  'RdYlBu',
#[Out]#  'RdYlGn',
#[Out]#  'Reds',
#[Out]#  'Spectral',
#[Out]#  'Wistia',
#[Out]#  'YlGn',
#[Out]#  'YlGnBu',
#[Out]#  'YlOrBr',
#[Out]#  'YlOrRd',
#[Out]#  'afmhot',
#[Out]#  'autumn',
#[Out]#  'binary',
#[Out]#  'bone',
#[Out]#  'brg',
#[Out]#  'bwr',
#[Out]#  'cool',
#[Out]#  'coolwarm',
#[Out]#  'copper',
#[Out]#  'cubehelix',
#[Out]#  'flag',
#[Out]#  'gist_earth',
#[Out]#  'gist_gray',
#[Out]#  'gist_heat',
#[Out]#  'gist_ncar',
#[Out]#  'gist_rainbow',
#[Out]#  'gist_stern',
#[Out]#  'gist_yarg',
#[Out]#  'gnuplot',
#[Out]#  'gnuplot2',
#[Out]#  'gray',
#[Out]#  'hot',
#[Out]#  'hsv',
#[Out]#  'jet',
#[Out]#  'nipy_spectral',
#[Out]#  'ocean',
#[Out]#  'pink',
#[Out]#  'prism',
#[Out]#  'rainbow',
#[Out]#  'seismic',
#[Out]#  'spring',
#[Out]#  'summer',
#[Out]#  'terrain',
#[Out]#  'winter',
#[Out]#  'Accent',
#[Out]#  'Dark2',
#[Out]#  'Paired',
#[Out]#  'Pastel1',
#[Out]#  'Pastel2',
#[Out]#  'Set1',
#[Out]#  'Set2',
#[Out]#  'Set3',
#[Out]#  'tab10',
#[Out]#  'tab20',
#[Out]#  'tab20b',
#[Out]#  'tab20c',
#[Out]#  'magma_r',
#[Out]#  'inferno_r',
#[Out]#  'plasma_r',
#[Out]#  'viridis_r',
#[Out]#  'cividis_r',
#[Out]#  'twilight_r',
#[Out]#  'twilight_shifted_r',
#[Out]#  'turbo_r',
#[Out]#  'Blues_r',
#[Out]#  'BrBG_r',
#[Out]#  'BuGn_r',
#[Out]#  'BuPu_r',
#[Out]#  'CMRmap_r',
#[Out]#  'GnBu_r',
#[Out]#  'Greens_r',
#[Out]#  'Greys_r',
#[Out]#  'OrRd_r',
#[Out]#  'Oranges_r',
#[Out]#  'PRGn_r',
#[Out]#  'PiYG_r',
#[Out]#  'PuBu_r',
#[Out]#  'PuBuGn_r',
#[Out]#  'PuOr_r',
#[Out]#  'PuRd_r',
#[Out]#  'Purples_r',
#[Out]#  'RdBu_r',
#[Out]#  'RdGy_r',
#[Out]#  'RdPu_r',
#[Out]#  'RdYlBu_r',
#[Out]#  'RdYlGn_r',
#[Out]#  'Reds_r',
#[Out]#  'Spectral_r',
#[Out]#  'Wistia_r',
#[Out]#  'YlGn_r',
#[Out]#  'YlGnBu_r',
#[Out]#  'YlOrBr_r',
#[Out]#  'YlOrRd_r',
#[Out]#  'afmhot_r',
#[Out]#  'autumn_r',
#[Out]#  'binary_r',
#[Out]#  'bone_r',
#[Out]#  'brg_r',
#[Out]#  'bwr_r',
#[Out]#  'cool_r',
#[Out]#  'coolwarm_r',
#[Out]#  'copper_r',
#[Out]#  'cubehelix_r',
#[Out]#  'flag_r',
#[Out]#  'gist_earth_r',
#[Out]#  'gist_gray_r',
#[Out]#  'gist_heat_r',
#[Out]#  'gist_ncar_r',
#[Out]#  'gist_rainbow_r',
#[Out]#  'gist_stern_r',
#[Out]#  'gist_yarg_r',
#[Out]#  'gnuplot_r',
#[Out]#  'gnuplot2_r',
#[Out]#  'gray_r',
#[Out]#  'hot_r',
#[Out]#  'hsv_r',
#[Out]#  'jet_r',
#[Out]#  'nipy_spectral_r',
#[Out]#  'ocean_r',
#[Out]#  'pink_r',
#[Out]#  'prism_r',
#[Out]#  'rainbow_r',
#[Out]#  'seismic_r',
#[Out]#  'spring_r',
#[Out]#  'summer_r',
#[Out]#  'terrain_r',
#[Out]#  'winter_r',
#[Out]#  'Accent_r',
#[Out]#  'Dark2_r',
#[Out]#  'Paired_r',
#[Out]#  'Pastel1_r',
#[Out]#  'Pastel2_r',
#[Out]#  'Set1_r',
#[Out]#  'Set2_r',
#[Out]#  'Set3_r',
#[Out]#  'tab10_r',
#[Out]#  'tab20_r',
#[Out]#  'tab20b_r',
#[Out]#  'tab20c_r']
pl.colormaps('viridis')
pl.colormaps()
#[Out]# ['magma',
#[Out]#  'inferno',
#[Out]#  'plasma',
#[Out]#  'viridis',
#[Out]#  'cividis',
#[Out]#  'twilight',
#[Out]#  'twilight_shifted',
#[Out]#  'turbo',
#[Out]#  'Blues',
#[Out]#  'BrBG',
#[Out]#  'BuGn',
#[Out]#  'BuPu',
#[Out]#  'CMRmap',
#[Out]#  'GnBu',
#[Out]#  'Greens',
#[Out]#  'Greys',
#[Out]#  'OrRd',
#[Out]#  'Oranges',
#[Out]#  'PRGn',
#[Out]#  'PiYG',
#[Out]#  'PuBu',
#[Out]#  'PuBuGn',
#[Out]#  'PuOr',
#[Out]#  'PuRd',
#[Out]#  'Purples',
#[Out]#  'RdBu',
#[Out]#  'RdGy',
#[Out]#  'RdPu',
#[Out]#  'RdYlBu',
#[Out]#  'RdYlGn',
#[Out]#  'Reds',
#[Out]#  'Spectral',
#[Out]#  'Wistia',
#[Out]#  'YlGn',
#[Out]#  'YlGnBu',
#[Out]#  'YlOrBr',
#[Out]#  'YlOrRd',
#[Out]#  'afmhot',
#[Out]#  'autumn',
#[Out]#  'binary',
#[Out]#  'bone',
#[Out]#  'brg',
#[Out]#  'bwr',
#[Out]#  'cool',
#[Out]#  'coolwarm',
#[Out]#  'copper',
#[Out]#  'cubehelix',
#[Out]#  'flag',
#[Out]#  'gist_earth',
#[Out]#  'gist_gray',
#[Out]#  'gist_heat',
#[Out]#  'gist_ncar',
#[Out]#  'gist_rainbow',
#[Out]#  'gist_stern',
#[Out]#  'gist_yarg',
#[Out]#  'gnuplot',
#[Out]#  'gnuplot2',
#[Out]#  'gray',
#[Out]#  'hot',
#[Out]#  'hsv',
#[Out]#  'jet',
#[Out]#  'nipy_spectral',
#[Out]#  'ocean',
#[Out]#  'pink',
#[Out]#  'prism',
#[Out]#  'rainbow',
#[Out]#  'seismic',
#[Out]#  'spring',
#[Out]#  'summer',
#[Out]#  'terrain',
#[Out]#  'winter',
#[Out]#  'Accent',
#[Out]#  'Dark2',
#[Out]#  'Paired',
#[Out]#  'Pastel1',
#[Out]#  'Pastel2',
#[Out]#  'Set1',
#[Out]#  'Set2',
#[Out]#  'Set3',
#[Out]#  'tab10',
#[Out]#  'tab20',
#[Out]#  'tab20b',
#[Out]#  'tab20c',
#[Out]#  'magma_r',
#[Out]#  'inferno_r',
#[Out]#  'plasma_r',
#[Out]#  'viridis_r',
#[Out]#  'cividis_r',
#[Out]#  'twilight_r',
#[Out]#  'twilight_shifted_r',
#[Out]#  'turbo_r',
#[Out]#  'Blues_r',
#[Out]#  'BrBG_r',
#[Out]#  'BuGn_r',
#[Out]#  'BuPu_r',
#[Out]#  'CMRmap_r',
#[Out]#  'GnBu_r',
#[Out]#  'Greens_r',
#[Out]#  'Greys_r',
#[Out]#  'OrRd_r',
#[Out]#  'Oranges_r',
#[Out]#  'PRGn_r',
#[Out]#  'PiYG_r',
#[Out]#  'PuBu_r',
#[Out]#  'PuBuGn_r',
#[Out]#  'PuOr_r',
#[Out]#  'PuRd_r',
#[Out]#  'Purples_r',
#[Out]#  'RdBu_r',
#[Out]#  'RdGy_r',
#[Out]#  'RdPu_r',
#[Out]#  'RdYlBu_r',
#[Out]#  'RdYlGn_r',
#[Out]#  'Reds_r',
#[Out]#  'Spectral_r',
#[Out]#  'Wistia_r',
#[Out]#  'YlGn_r',
#[Out]#  'YlGnBu_r',
#[Out]#  'YlOrBr_r',
#[Out]#  'YlOrRd_r',
#[Out]#  'afmhot_r',
#[Out]#  'autumn_r',
#[Out]#  'binary_r',
#[Out]#  'bone_r',
#[Out]#  'brg_r',
#[Out]#  'bwr_r',
#[Out]#  'cool_r',
#[Out]#  'coolwarm_r',
#[Out]#  'copper_r',
#[Out]#  'cubehelix_r',
#[Out]#  'flag_r',
#[Out]#  'gist_earth_r',
#[Out]#  'gist_gray_r',
#[Out]#  'gist_heat_r',
#[Out]#  'gist_ncar_r',
#[Out]#  'gist_rainbow_r',
#[Out]#  'gist_stern_r',
#[Out]#  'gist_yarg_r',
#[Out]#  'gnuplot_r',
#[Out]#  'gnuplot2_r',
#[Out]#  'gray_r',
#[Out]#  'hot_r',
#[Out]#  'hsv_r',
#[Out]#  'jet_r',
#[Out]#  'nipy_spectral_r',
#[Out]#  'ocean_r',
#[Out]#  'pink_r',
#[Out]#  'prism_r',
#[Out]#  'rainbow_r',
#[Out]#  'seismic_r',
#[Out]#  'spring_r',
#[Out]#  'summer_r',
#[Out]#  'terrain_r',
#[Out]#  'winter_r',
#[Out]#  'Accent_r',
#[Out]#  'Dark2_r',
#[Out]#  'Paired_r',
#[Out]#  'Pastel1_r',
#[Out]#  'Pastel2_r',
#[Out]#  'Set1_r',
#[Out]#  'Set2_r',
#[Out]#  'Set3_r',
#[Out]#  'tab10_r',
#[Out]#  'tab20_r',
#[Out]#  'tab20b_r',
#[Out]#  'tab20c_r']
pl.matplotlib.cm.viridis
#[Out]# <matplotlib.colors.ListedColormap at 0x2baa9329d940>
pl.matplotlib.cm.viridis.copy()
#[Out]# <matplotlib.colors.ListedColormap at 0x2bab379a7ee0>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='asinh', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm, cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6e3402b0>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='log', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', norm=norm, cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6eff5970>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6f1d3310>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6f394f10>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6f651fd0>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)
ax2.xaxis.set_ticks([])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6f7d2730>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)
ax2.xaxis.set_tickabels([])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)
ax2.xaxis.set_ticklabels([])

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6fc58fd0>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (10, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2bab6fe34e50>
mask = signal_mask.sum(axis=0)
collapse1 = signal_mask.sum(axis = 1)
collapse2 = signal_mask.sum(axis = 2)

vmax = np.max([mask.max(), collapse1.max(), collapse2.max()])
norm = simple_norm(mask, stretch='linear', max_cut=vmax, min_cut=0.1)

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w')

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize = (9, 10), constrained_layout=True)

ax1 = plt.subplot(323, aspect = 1)
im1 = ax1.imshow(mask, origin='lower', cmap=cm, vmin=0.1) # cmap='inferno', aspect='auto'
ax1.tick_params(direction='in', color='#FFFFFF')
pl.colorbar(mappable=im1)

ax2 = plt.subplot(321, sharex = ax1, adjustable='box')
im2 = ax2.imshow(collapse1, origin='lower', norm=norm, cmap=cm) # cmap='gray'
ax2.set_ylabel('Velocity (km/s)')
ax2.tick_params(direction='in', color='#FFFFFF')
ax2.set_aspect(collapse1.shape[1]/collapse1.shape[0])
pl.colorbar(mappable=im2)

ax3 = plt.subplot(324, sharey = ax1)
im3 = ax3.imshow(collapse2.T, origin='lower', norm=norm, cmap=cm) # cmap='gray', aspect='auto'
ax3.set_xlabel('Velocity (km/s)')
ax3.tick_params(direction='in', color='#FFFFFF')
ax3.set_aspect(collapse2.shape[0]/collapse2.shape[1])

#pl.colorbar(cax=[0.95,0.5,0.1,0.8], mappable=im3)
pl.colorbar(mappable=im3)

# fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
#                         gridspec_kw=dict(height_ratios=[1, 3],
#                                          width_ratios=[3, 1]))
# axs[0, 1].set_visible(False)
# axs[0, 0].set_box_aspect(1/3)
# axs[1, 0].set_box_aspect(1)
# axs[1, 1].set_box_aspect(3/1)

# x, y = np.random.randn(2, 400) * [[.5], [180]]
# axs[1, 0].imshow(mask) # origin='lower', adjustable='box'
# axs[0, 0].imshow(collapse1)
# axs[1, 1].imshow(collapse2.T)

# plt.show()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2babda5a90d0>
