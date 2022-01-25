########################################################
# Started Logging At: 2022-01-19 16:24:17
########################################################
########################################################
# # Started Logging At: 2022-01-19 16:24:18
########################################################
# Make thioformaldehyde (H2CS) cubes for inspection with DS9
from astroquery.splatalogue import Splatalogue
from spectral_cube import SpectralCube
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as nd
########################################################
# Started Logging At: 2022-01-19 16:24:29
########################################################
########################################################
# # Started Logging At: 2022-01-19 16:24:29
########################################################
# Make thioformaldehyde (H2CS) cubes for inspection with DS9
from astroquery.splatalogue import Splatalogue
from spectral_cube import SpectralCube
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as nd
tbl = Splatalogue.query_lines(85*u.GHz, 155*u.GHz, chemical_name=' H2CS ', line_lists=['JPL'],
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
tbl = Splatalogue.query_lines(239.9*u.GHz, 275*u.GHz, chemical_name=' H2CS ', line_lists=['JPL'], 
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '137_spw69'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=137.371051*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_404-303.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '102_spw29'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=103.04022*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_303-202.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '102_spw23'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=101.47762*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_313-212.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '103_spw31'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=104.61704*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_312-211.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '135_spw47'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=135.297811*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_414-313.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '139_spw71'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=139.48341*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_413-312.fits')
tbl = Splatalogue.query_lines(80*u.GHz, 90*u.GHz, chemical_name=' NH2D ', line_lists=['JPL'], 
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '87_spw25'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
nh2dcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=85.9247829*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# nh2dcube.write('NH2D_hyperfine.fits')
def get_noise_map(cube_noise):
    cube_sclip = cube_noise.sigma_clip_spectrally(3) # Clip values above 3-sigma 
    mad_std_spectrum_sclip = cube_sclip.mad_std(axis=(1, 2))
    plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
             drawstyle='steps-mid')
    plt.xlabel('Velocity (km/s)')
    plt.ylabel(r' Noise standard deviation $\sigma$ (K)')
#     plt.ylim([0., 0.30]) # Best to extend the range to 0.
#     plt.axhline(0.25, linestyle='--', color='k', linewidth=3, label='A priori noise expectation')
    plt.legend(frameon=True)
    cube_sclip_slice = cube_sclip.spectral_slab(50*u.km/u.s, 80*u.km/u.s)
    mad_std_map_sclip = cube_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
    return mad_std_map_sclip

def get_signal_mask_scipy(cube_signal, mad_std_map_sclip):
    '''Please pass an already-masked cube to cube_signal'''
    # Make a low and high mask
    low_snr_mask = (cube_signal > 3 * mad_std_map_sclip).include()
    high_snr_mask = (cube_signal > 10 * mad_std_map_sclip).include()
    low_snr_mask = low_snr_mask.compute() # We need to convert from a dask array to a numpy array.
    high_snr_mask = high_snr_mask.compute()
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
    return signal_mask
# Get overall cube
cube = SpectralCube.read('H2CS_414-313.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(85, 512, 512) and unit=Jy / beam and chunk size (85, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     85  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      80.460 km / s
subcube_low = cube.spectral_slab(-10*u.km/u.s, 25*u.km/u.s)
subcube_high = cube.spectral_slab(25*u.km/u.s, 50*u.km/u.s)
# Get noise map
mad_std_map_sclip_low = get_noise_map(subcube_low)
mad_std_map_sclip_high = get_noise_map(subcube_high)
plain_mask_low = subcube_low >= 3 * mad_std_map_sclip_low
plain_masked_slab_low = subcube_low.with_mask(plain_mask_low)
plain_mask_high = subcube_high >= 3 * mad_std_map_sclip_high
plain_masked_slab_high = subcube_high.with_mask(plain_mask_high)
signal_mask_low = get_signal_mask_scipy(plain_masked_slab_low, mad_std_map_sclip_low)
masked_cube_low = plain_masked_slab_low.with_mask(signal_mask_low)
signal_mask_high = get_signal_mask_scipy(plain_masked_slab_high, mad_std_map_sclip_high)
masked_cube_high = plain_masked_slab_high.with_mask(signal_mask_high)
# masked_moment0 = masked_cube.moment0() # masked_cubes[n]

# ax = plt.subplot(projection=masked_moment0.wcs)
# im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
# cbar = plt.colorbar(im)
# cbar.set_label('Integrated Intensity (K km/s)')
# ax.set_ylabel('Declination')
# ax.set_xlabel('Right Ascension')
# masked_moment1 = masked_cube.moment1() # masked_cubes[n]

# ax = plt.subplot(projection=masked_moment1.wcs)
# im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
# cbar = plt.colorbar(im)
# cbar.set_label('Centroid (km/s)')

# ax.set_ylabel('Declination')
# ax.set_xlabel('Right Ascension')
# from astroquery.splatalogue import Splatalogue
# from astropy import units as u
# from spectral_cube import SpectralCube
# results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# tbl = Splatalogue.query_lines(109*u.GHz, 111*u.GHz, chemical_name=' NH2D ', line_lists=['JPL'],
#                               energy_max=150, energy_type='eu_k', show_qn_code=True)
# tbl
# freq_spw = '110_spw29'
# fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
# cube = SpectralCube.read(fn, format='casa_image')
# cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
#                                    rest_value=110.1536299*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# # cube.write('test_nh2d.fits')
# cube = SpectralCube.read('test_nh2d.fits')
# cube.max(axis=0).quicklook()
# # plt.plot(cube.max(axis=0)) # this is not right but it does look cool
# mxspec = cube.max(axis=(1,2))
# import pylab as pl
# import pyspeckit
# kspectrum_ps = pyspeckit.Spectrum(xarr = mxspec.spectral_axis, data = mxspec.value)
# Let's get masked moment maps for each of the cubes
masked_moment0_low = masked_cube_low.moment0()
masked_moment0_high = masked_cube_high.moment0()
mom0s = [masked_moment0_low, masked_moment0_high]
masked_cube_low.with_spectral_unit(u.km/u.s)
#[Out]# DaskVaryingResolutionSpectralCube with shape=(34, 512, 512) and unit=Jy / beam and chunk size (34, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     34  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      25.280 km / s
collapse1_low = masked_cube_low.with_spectral_unit(u.km/u.s).mean(axis = 1) # Collapse along one spatial axis
collapse2_low = masked_cube_low.with_spectral_unit(u.km/u.s).mean(axis = 2) # Collapse along the other spatial axis
collapse1_high = masked_cube_high.with_spectral_unit(u.km/u.s).mean(axis = 1) # Collapse along one spatial axis
collapse2_high = masked_cube_high.with_spectral_unit(u.km/u.s).mean(axis = 2) # Collapse along the other spatial axis
collapse1s = [collapse1_low, collapse1_high]
collapse2s = [collapse2_low, collapse2_high]
# collapse1.wcs
from mpl_toolkits.axes_grid1 import make_axes_locatable

i = 0
for i in range(len(mom0s)):
    masked_moment0 = mom0s[i]
    collapse1 = collapse1s[i]
    collapse2 = collapse2s[i]
    if i == 0:
        plot_type = 'Low'
    elif i == 1:
        plot_type = 'High'
    
    # Initalize figure
    fig = plt.figure(figsize = (10, 10), constrained_layout=True)
    
    # Plot moment 0 map
    ax1 = plt.subplot(323, projection = masked_moment0.wcs)
    ax1.set_aspect(1.)
    im1 = ax1.imshow(masked_moment0.value, origin='lower', cmap='inferno', aspect='auto')
    ax1.set_ylabel('Declination')
    ax1.set_xlabel('Right Ascension')
#     plt.title(f'H2CS $4(1, 4)-3(1, 3)$ Moment 0 Map ({plot_type} Velocity)')
    ax1.tick_params(direction='in', color='#FFFFFF')
    ax1.set_facecolor('#000000')
    
#     cax = plt.subplot(325)
#     cbar = plt.colorbar(im1, cax = cax, orientation = 'horizontal')
#     cbar.set_label('Integrated Intensity (K km/s)')

    # Plot PV projection along one spatial axis (TOP LEFT)
#     ax2 = divider.append_axes("top", size="5%", axes_class = ax1.Axes)
    ax2 = plt.subplot(321, projection = collapse1.wcs, sharex = ax1)
#     ax2 = plt.subplot(321, sharex = ax1)
    im2 = ax2.imshow(collapse1.value, origin='lower', cmap='gray', aspect='auto')
    ax2.coords[1].set_format_unit(u.km / u.s)
    ax2.set_ylabel('Velocity (km/s)')
    ax2.coords[0].set_ticklabel_visible(False)
    ax2.tick_params(direction='in', color='#FFFFFF')
    ax2.set_facecolor('#000000')

    # Plot PV projection along the other spatial axis (BOTTOM RIGHT)
#     ax3 = divider.append_axes("right", size="5%", axes_class = ax1.Axes)
    ax3 = plt.subplot(324, projection = collapse2.wcs.sub([2,1]), sharey = ax1)
    im3 = ax3.imshow(collapse2.T.value, origin='lower', cmap='gray', aspect='auto')
    ax3.coords[0].set_format_unit(u.km / u.s)
    ax3.set_xlabel('Velocity (km/s)')
    ax3.coords[1].set_ticklabel_visible(False)
    ax3.tick_params(direction='in', color='#FFFFFF')
    ax3.set_facecolor('#000000')
    
    # Trying to remove axis labels
#     plt.setp(ax3.get_yticklabels(), visible=False) # This didn't remove anything
#     ax3.get_yaxis().set_visible(False) # This also didn't remove anything
#     ax3.yaxis.set_visible(False) # Also does nothing
#     plt.axis('off') # Not what I want, but works

    # Checking axes
#     print(f"Moment map x limits: {ax1.get_xlim()}")
#     print(f"Upper PV projection x limits: {ax2.get_xlim()}")
#     print(f"Right PV projection x limits: {ax3.get_xlim()}")
#     print(f"Moment map y limits: {ax1.get_ylim()}")
#     print(f"Upper PV projection y limits: {ax2.get_ylim()}")
#     print(f"Right PV projection y limits: {ax3.get_ylim()}")

    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/thioformaldehyde/prelim_mom0_PV_projection.pdf', bbox_inches='tight')
    plt.show()
    i += 1
# # STILL NEED TO:
# # Match coordinate axes
# # Switch axis numbers and labels for PV over RA
# # Make plot labels nice

# i = 0
# for i in range(len(mom0s)):
#     masked_moment0 = mom0s[i]
#     collapse1 = collapse1s[i]
#     collapse2 = collapse2s[i]
#     if i == 0:
#         plot_type = 'Low'
#     elif i == 1:
#         plot_type = 'High'
    
#     # Initalize figure
# #     plt.subplots(2, 2, figsize = (10, 10), sharex='col', sharey='row')
#     fig = plt.figure(figsize = (10, 10)) # tight_layout = True
    
#     # Plot moment 0 map
#     ax1 = fig.add_subplot(2, 2, 3, projection = masked_moment0.wcs)
#     im1 = plt.imshow(masked_moment0.value, origin='lower',  cmap='inferno')
#     ax1.set_ylabel('Declination')
#     ax1.set_xlabel('Right Ascension')
#     plt.title(f'H2CS $4(1, 4)-3(1, 3)$ Moment 0 Map ({plot_type} Velocity)')

#     # Plot PV projection along one spatial axis
#     ax2 = fig.add_subplot(2, 2, 1, projection = collapse1.wcs, sharex = ax1)
#     im2 = plt.imshow(collapse1.value, origin='lower',  cmap='inferno', aspect='auto') # aspect='auto'
# #     ax2.set_ylabel('Velocity')
# #     ax2.set_xlabel('Right Ascension')
# #     ax2.set_xticks(ax1.get_xticks())
# #     ax2.set_xticklabels(ax1.get_xticklabels())
#     plt.title('Position-velocity projection along dec')

#     # Plot PV projection along the other spatial axis
#     ax3 = fig.add_subplot(2, 2, 4, projection = collapse2.T.wcs, sharey = ax1)
#     im3 = plt.imshow(collapse2.T.value, origin='lower', cmap='inferno', aspect='auto')
# #     ax3.set_ylabel('Velocity')
# #     ax3.set_xlabel('Declination')
#     plt.title('Position-velocity projection along RA')

#     cbar = plt.colorbar(im1)
#     cbar.set_label('Integrated Intensity (K km/s)')

#     # plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/lineID/plots/unlabeled/{freq_spw}.pdf')
#     plt.show()
#     i += 1
# Get overall cube
cube = SpectralCube.read('NH2D_hyperfine.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(27, 512, 512) and unit=Jy / beam and chunk size (27, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     27  type_s: VRAD      unit_s: km / s  range:       -9.493 km / s:      79.083 km / s
subcube_low = cube.spectral_slab(-10*u.km/u.s, 25*u.km/u.s)
subcube_high = cube.spectral_slab(25*u.km/u.s, 50*u.km/u.s)
# Get noise map
mad_std_map_sclip_low = get_noise_map(subcube_low)
mad_std_map_sclip_high = get_noise_map(subcube_high)
plain_mask_low = subcube_low >= 3 * mad_std_map_sclip_low
plain_masked_slab_low = subcube_low.with_mask(plain_mask_low)
plain_mask_high = subcube_high >= 3 * mad_std_map_sclip_high
plain_masked_slab_high = subcube_high.with_mask(plain_mask_high)
signal_mask_low = get_signal_mask_scipy(plain_masked_slab_low, mad_std_map_sclip_low)
masked_cube_low = plain_masked_slab_low.with_mask(signal_mask_low)
signal_mask_high = get_signal_mask_scipy(plain_masked_slab_high, mad_std_map_sclip_high)
masked_cube_high = plain_masked_slab_high.with_mask(signal_mask_high)
# Let's get masked moment maps for each of the cubes
masked_moment0_low = masked_cube_low.moment0()
masked_moment0_high = masked_cube_high.moment0()
mom0s = [masked_moment0_low, masked_moment0_high]
collapse1_low = masked_cube_low.mean(axis = 1) # Collapse along one spatial axis
collapse2_low = masked_cube_low.mean(axis = 2) # Collapse along the other spatial axis
collapse1_high = masked_cube_high.mean(axis = 1) # Collapse along one spatial axis
collapse2_high = masked_cube_high.mean(axis = 2) # Collapse along the other spatial axis
collapse1s = [collapse1_low, collapse1_high]
collapse2s = [collapse2_low, collapse2_high]
i = 0
for i in range(len(mom0s)):
    masked_moment0 = mom0s[i]
    collapse1 = collapse1s[i]
    collapse2 = collapse2s[i]
    if i == 0:
        plot_type = 'Low'
    elif i == 1:
        plot_type = 'High'
    
    # Initalize figure
    fig = plt.figure(1, figsize = (10, 10), tight_layout = False)

    # Plot moment 0 map
    ax = plt.subplot(223, projection = masked_moment0.wcs)
    im1 = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
    ax.set_ylabel('Declination')
    ax.set_xlabel('Right Ascension')
    plt.title(f'NH2D 85 GHz Hyperfine Moment 0 Map ({plot_type} Velocity)')

    # Plot PV projection along one spatial axis
    ax = plt.subplot(221, projection = collapse1.wcs)
    im2 = ax.imshow(collapse1.value, origin='lower', cmap='inferno', aspect='auto')
    ax.set_ylabel('Velocity')
    ax.set_xlabel('Right Ascension')
    plt.title('Position-velocity projection 1')

    # Plot PV projection along the other spatial axis
    ax = plt.subplot(224, projection = collapse2.T.wcs)
    im3 = ax.imshow(collapse2.T.value, origin='lower', cmap='inferno', aspect='auto')
    ax.set_ylabel('Velocity')
    ax.set_xlabel('Declination')
    plt.title('Position-velocity projection 2')

    cbar = plt.colorbar(im1)
    cbar.set_label('Integrated Intensity (K km/s)')

    # plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/lineID/plots/unlabeled/{freq_spw}.pdf')
    plt.show()
    i += 1
# Make thioformaldehyde (H2CS) cubes for inspection with DS9
from astroquery.splatalogue import Splatalogue
from spectral_cube import SpectralCube
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as nd
plt.rcParams['figure.facecolor'] = 'w'
tbl = Splatalogue.query_lines(85*u.GHz, 155*u.GHz, chemical_name=' H2CS ', line_lists=['JPL'],
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
tbl = Splatalogue.query_lines(239.9*u.GHz, 275*u.GHz, chemical_name=' H2CS ', line_lists=['JPL'], 
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '137_spw69'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=137.371051*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_404-303.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '102_spw29'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=103.04022*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_303-202.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '102_spw23'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=101.47762*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_313-212.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '103_spw31'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=104.61704*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_312-211.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '135_spw47'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=135.297811*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_414-313.fits')
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '139_spw71'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
h2cscube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=139.48341*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# h2cscube.write('H2CS_413-312.fits')
tbl = Splatalogue.query_lines(80*u.GHz, 90*u.GHz, chemical_name=' NH2D ', line_lists=['JPL'], 
                              energy_max=150, energy_type='eu_k', show_qn_code=True)
tbl.show_in_notebook()
#[Out]# <IPython.core.display.HTML object>
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '87_spw25'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube = SpectralCube.read(fn, format='casa_image')
nh2dcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
                                   rest_value=85.9247829*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# nh2dcube.write('NH2D_hyperfine.fits')
def get_noise_map(cube_noise):
    cube_sclip = cube_noise.sigma_clip_spectrally(3) # Clip values above 3-sigma 
    mad_std_spectrum_sclip = cube_sclip.mad_std(axis=(1, 2))
    plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, 
             drawstyle='steps-mid')
    plt.xlabel('Velocity (km/s)')
    plt.ylabel(r' Noise standard deviation $\sigma$ (K)')
#     plt.ylim([0., 0.30]) # Best to extend the range to 0.
#     plt.axhline(0.25, linestyle='--', color='k', linewidth=3, label='A priori noise expectation')
    plt.legend(frameon=True)
    cube_sclip_slice = cube_sclip.spectral_slab(50*u.km/u.s, 80*u.km/u.s)
    mad_std_map_sclip = cube_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
    return mad_std_map_sclip

def get_signal_mask_scipy(cube_signal, mad_std_map_sclip):
    '''Please pass an already-masked cube to cube_signal'''
    # Make a low and high mask
    low_snr_mask = (cube_signal > 3 * mad_std_map_sclip).include()
    high_snr_mask = (cube_signal > 10 * mad_std_map_sclip).include()
    low_snr_mask = low_snr_mask.compute() # We need to convert from a dask array to a numpy array.
    high_snr_mask = high_snr_mask.compute()
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
    return signal_mask
# Get overall cube
cube = SpectralCube.read('H2CS_414-313.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(85, 512, 512) and unit=Jy / beam and chunk size (85, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     85  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      80.460 km / s
subcube_low = cube.spectral_slab(-10*u.km/u.s, 25*u.km/u.s)
subcube_high = cube.spectral_slab(25*u.km/u.s, 50*u.km/u.s)
# Get noise map
mad_std_map_sclip_low = get_noise_map(subcube_low)
mad_std_map_sclip_high = get_noise_map(subcube_high)
plain_mask_low = subcube_low >= 3 * mad_std_map_sclip_low
plain_masked_slab_low = subcube_low.with_mask(plain_mask_low)
plain_mask_high = subcube_high >= 3 * mad_std_map_sclip_high
plain_masked_slab_high = subcube_high.with_mask(plain_mask_high)
signal_mask_low = get_signal_mask_scipy(plain_masked_slab_low, mad_std_map_sclip_low)
masked_cube_low = plain_masked_slab_low.with_mask(signal_mask_low)
signal_mask_high = get_signal_mask_scipy(plain_masked_slab_high, mad_std_map_sclip_high)
masked_cube_high = plain_masked_slab_high.with_mask(signal_mask_high)
# masked_moment0 = masked_cube.moment0() # masked_cubes[n]

# ax = plt.subplot(projection=masked_moment0.wcs)
# im = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
# cbar = plt.colorbar(im)
# cbar.set_label('Integrated Intensity (K km/s)')
# ax.set_ylabel('Declination')
# ax.set_xlabel('Right Ascension')
# masked_moment1 = masked_cube.moment1() # masked_cubes[n]

# ax = plt.subplot(projection=masked_moment1.wcs)
# im = ax.imshow(masked_moment1.value, origin='lower', cmap='coolwarm')
# cbar = plt.colorbar(im)
# cbar.set_label('Centroid (km/s)')

# ax.set_ylabel('Declination')
# ax.set_xlabel('Right Ascension')
# from astroquery.splatalogue import Splatalogue
# from astropy import units as u
# from spectral_cube import SpectralCube
# results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# tbl = Splatalogue.query_lines(109*u.GHz, 111*u.GHz, chemical_name=' NH2D ', line_lists=['JPL'],
#                               energy_max=150, energy_type='eu_k', show_qn_code=True)
# tbl
# freq_spw = '110_spw29'
# fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
# cube = SpectralCube.read(fn, format='casa_image')
# cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', 
#                                    rest_value=110.1536299*u.GHz).spectral_slab(-10*u.km/u.s, 80*u.km/u.s)
# # cube.write('test_nh2d.fits')
# cube = SpectralCube.read('test_nh2d.fits')
# cube.max(axis=0).quicklook()
# # plt.plot(cube.max(axis=0)) # this is not right but it does look cool
# mxspec = cube.max(axis=(1,2))
# import pylab as pl
# import pyspeckit
# kspectrum_ps = pyspeckit.Spectrum(xarr = mxspec.spectral_axis, data = mxspec.value)
# Let's get masked moment maps for each of the cubes
masked_moment0_low = masked_cube_low.moment0()
masked_moment0_high = masked_cube_high.moment0()
mom0s = [masked_moment0_low, masked_moment0_high]
masked_cube_low.with_spectral_unit(u.km/u.s)
#[Out]# DaskVaryingResolutionSpectralCube with shape=(34, 512, 512) and unit=Jy / beam and chunk size (34, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     34  type_s: VRAD      unit_s: km / s  range:      -10.425 km / s:      25.280 km / s
collapse1_low = masked_cube_low.with_spectral_unit(u.km/u.s).mean(axis = 1) # Collapse along one spatial axis
collapse2_low = masked_cube_low.with_spectral_unit(u.km/u.s).mean(axis = 2) # Collapse along the other spatial axis
collapse1_high = masked_cube_high.with_spectral_unit(u.km/u.s).mean(axis = 1) # Collapse along one spatial axis
collapse2_high = masked_cube_high.with_spectral_unit(u.km/u.s).mean(axis = 2) # Collapse along the other spatial axis
collapse1s = [collapse1_low, collapse1_high]
collapse2s = [collapse2_low, collapse2_high]
# collapse1.wcs
from mpl_toolkits.axes_grid1 import make_axes_locatable

i = 0
for i in range(len(mom0s)):
    masked_moment0 = mom0s[i]
    collapse1 = collapse1s[i]
    collapse2 = collapse2s[i]
    if i == 0:
        plot_type = 'Low'
    elif i == 1:
        plot_type = 'High'
    
    # Initalize figure
    fig = plt.figure(figsize = (10, 10), constrained_layout=True)
    
    # Plot moment 0 map
    ax1 = plt.subplot(323, projection = masked_moment0.wcs)
    ax1.set_aspect(1.)
    im1 = ax1.imshow(masked_moment0.value, origin='lower', cmap='inferno', aspect='auto')
    ax1.set_ylabel('Declination')
    ax1.set_xlabel('Right Ascension')
#     plt.title(f'H2CS $4(1, 4)-3(1, 3)$ Moment 0 Map ({plot_type} Velocity)')
    ax1.tick_params(direction='in', color='#FFFFFF')
    ax1.set_facecolor('#000000')
    
#     cax = plt.subplot(325)
#     cbar = plt.colorbar(im1, cax = cax, orientation = 'horizontal')
#     cbar.set_label('Integrated Intensity (K km/s)')

    # Plot PV projection along one spatial axis (TOP LEFT)
#     ax2 = divider.append_axes("top", size="5%", axes_class = ax1.Axes)
    ax2 = plt.subplot(321, projection = collapse1.wcs, sharex = ax1)
#     ax2 = plt.subplot(321, sharex = ax1)
    im2 = ax2.imshow(collapse1.value, origin='lower', cmap='gray', aspect='auto')
    ax2.coords[1].set_format_unit(u.km / u.s)
    ax2.set_ylabel('Velocity (km/s)')
    ax2.coords[0].set_ticklabel_visible(False)
    ax2.tick_params(direction='in', color='#FFFFFF')
    ax2.set_facecolor('#000000')

    # Plot PV projection along the other spatial axis (BOTTOM RIGHT)
#     ax3 = divider.append_axes("right", size="5%", axes_class = ax1.Axes)
    ax3 = plt.subplot(324, projection = collapse2.wcs.sub([2,1]), sharey = ax1)
    im3 = ax3.imshow(collapse2.T.value, origin='lower', cmap='gray', aspect='auto')
    ax3.coords[0].set_format_unit(u.km / u.s)
    ax3.set_xlabel('Velocity (km/s)')
    ax3.coords[1].set_ticklabel_visible(False)
    ax3.tick_params(direction='in', color='#FFFFFF')
    ax3.set_facecolor('#000000')
    
    # Trying to remove axis labels
#     plt.setp(ax3.get_yticklabels(), visible=False) # This didn't remove anything
#     ax3.get_yaxis().set_visible(False) # This also didn't remove anything
#     ax3.yaxis.set_visible(False) # Also does nothing
#     plt.axis('off') # Not what I want, but works

    # Checking axes
#     print(f"Moment map x limits: {ax1.get_xlim()}")
#     print(f"Upper PV projection x limits: {ax2.get_xlim()}")
#     print(f"Right PV projection x limits: {ax3.get_xlim()}")
#     print(f"Moment map y limits: {ax1.get_ylim()}")
#     print(f"Upper PV projection y limits: {ax2.get_ylim()}")
#     print(f"Right PV projection y limits: {ax3.get_ylim()}")

    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/thioformaldehyde/prelim_mom0_PV_projection.pdf', bbox_inches='tight')
    plt.show()
    i += 1
# # STILL NEED TO:
# # Match coordinate axes
# # Switch axis numbers and labels for PV over RA
# # Make plot labels nice

# i = 0
# for i in range(len(mom0s)):
#     masked_moment0 = mom0s[i]
#     collapse1 = collapse1s[i]
#     collapse2 = collapse2s[i]
#     if i == 0:
#         plot_type = 'Low'
#     elif i == 1:
#         plot_type = 'High'
    
#     # Initalize figure
# #     plt.subplots(2, 2, figsize = (10, 10), sharex='col', sharey='row')
#     fig = plt.figure(figsize = (10, 10)) # tight_layout = True
    
#     # Plot moment 0 map
#     ax1 = fig.add_subplot(2, 2, 3, projection = masked_moment0.wcs)
#     im1 = plt.imshow(masked_moment0.value, origin='lower',  cmap='inferno')
#     ax1.set_ylabel('Declination')
#     ax1.set_xlabel('Right Ascension')
#     plt.title(f'H2CS $4(1, 4)-3(1, 3)$ Moment 0 Map ({plot_type} Velocity)')

#     # Plot PV projection along one spatial axis
#     ax2 = fig.add_subplot(2, 2, 1, projection = collapse1.wcs, sharex = ax1)
#     im2 = plt.imshow(collapse1.value, origin='lower',  cmap='inferno', aspect='auto') # aspect='auto'
# #     ax2.set_ylabel('Velocity')
# #     ax2.set_xlabel('Right Ascension')
# #     ax2.set_xticks(ax1.get_xticks())
# #     ax2.set_xticklabels(ax1.get_xticklabels())
#     plt.title('Position-velocity projection along dec')

#     # Plot PV projection along the other spatial axis
#     ax3 = fig.add_subplot(2, 2, 4, projection = collapse2.T.wcs, sharey = ax1)
#     im3 = plt.imshow(collapse2.T.value, origin='lower', cmap='inferno', aspect='auto')
# #     ax3.set_ylabel('Velocity')
# #     ax3.set_xlabel('Declination')
#     plt.title('Position-velocity projection along RA')

#     cbar = plt.colorbar(im1)
#     cbar.set_label('Integrated Intensity (K km/s)')

#     # plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/lineID/plots/unlabeled/{freq_spw}.pdf')
#     plt.show()
#     i += 1
# Get overall cube
cube = SpectralCube.read('NH2D_hyperfine.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(27, 512, 512) and unit=Jy / beam and chunk size (27, 512, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     27  type_s: VRAD      unit_s: km / s  range:       -9.493 km / s:      79.083 km / s
subcube_low = cube.spectral_slab(-10*u.km/u.s, 25*u.km/u.s)
subcube_high = cube.spectral_slab(25*u.km/u.s, 50*u.km/u.s)
# Get noise map
mad_std_map_sclip_low = get_noise_map(subcube_low)
mad_std_map_sclip_high = get_noise_map(subcube_high)
plain_mask_low = subcube_low >= 3 * mad_std_map_sclip_low
plain_masked_slab_low = subcube_low.with_mask(plain_mask_low)
plain_mask_high = subcube_high >= 3 * mad_std_map_sclip_high
plain_masked_slab_high = subcube_high.with_mask(plain_mask_high)
signal_mask_low = get_signal_mask_scipy(plain_masked_slab_low, mad_std_map_sclip_low)
masked_cube_low = plain_masked_slab_low.with_mask(signal_mask_low)
signal_mask_high = get_signal_mask_scipy(plain_masked_slab_high, mad_std_map_sclip_high)
masked_cube_high = plain_masked_slab_high.with_mask(signal_mask_high)
# Let's get masked moment maps for each of the cubes
masked_moment0_low = masked_cube_low.moment0()
masked_moment0_high = masked_cube_high.moment0()
mom0s = [masked_moment0_low, masked_moment0_high]
collapse1_low = masked_cube_low.mean(axis = 1) # Collapse along one spatial axis
collapse2_low = masked_cube_low.mean(axis = 2) # Collapse along the other spatial axis
collapse1_high = masked_cube_high.mean(axis = 1) # Collapse along one spatial axis
collapse2_high = masked_cube_high.mean(axis = 2) # Collapse along the other spatial axis
collapse1s = [collapse1_low, collapse1_high]
collapse2s = [collapse2_low, collapse2_high]
i = 0
for i in range(len(mom0s)):
    masked_moment0 = mom0s[i]
    collapse1 = collapse1s[i]
    collapse2 = collapse2s[i]
    if i == 0:
        plot_type = 'Low'
    elif i == 1:
        plot_type = 'High'
    
    # Initalize figure
    fig = plt.figure(1, figsize = (10, 10), tight_layout = False)

    # Plot moment 0 map
    ax = plt.subplot(223, projection = masked_moment0.wcs)
    im1 = ax.imshow(masked_moment0.value, origin='lower', cmap='inferno')
    ax.set_ylabel('Declination')
    ax.set_xlabel('Right Ascension')
    plt.title(f'NH2D 85 GHz Hyperfine Moment 0 Map ({plot_type} Velocity)')

    # Plot PV projection along one spatial axis
    ax = plt.subplot(221, projection = collapse1.wcs)
    im2 = ax.imshow(collapse1.value, origin='lower', cmap='inferno', aspect='auto')
    ax.set_ylabel('Velocity')
    ax.set_xlabel('Right Ascension')
    plt.title('Position-velocity projection 1')

    # Plot PV projection along the other spatial axis
    ax = plt.subplot(224, projection = collapse2.T.wcs)
    im3 = ax.imshow(collapse2.T.value, origin='lower', cmap='inferno', aspect='auto')
    ax.set_ylabel('Velocity')
    ax.set_xlabel('Declination')
    plt.title('Position-velocity projection 2')

    cbar = plt.colorbar(im1)
    cbar.set_label('Integrated Intensity (K km/s)')

    # plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/lineID/plots/unlabeled/{freq_spw}.pdf')
    plt.show()
    i += 1
