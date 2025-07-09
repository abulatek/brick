########################################################
# Started Logging At: 2021-12-01 14:38:13
########################################################
########################################################
# # Started Logging At: 2021-12-01 14:38:17
########################################################
########################################################
# Started Logging At: 2021-12-01 14:38:17
########################################################
########################################################
# # Started Logging At: 2021-12-01 14:38:17
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
def retrieve_cube(results, freq_spw):
    '''Get methyl cyanide (target molecule) cube'''
    fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
    ch3cncube = SpectralCube.read(fn, format='casa_image')
    ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN',
                                                                                        fmin=fmin, 
                                                                                        fmax=fmax, 
                                                                                        catalog='JPL')
    # We're readying the partition function for use with temperature map later!
    ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
    ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data
    imshow(noise_map, origin='lower')
    return noise_map
noise_map = generate_noise_map()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        noise_map_int = noise_map*u.K*(channel_width)

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

        # Append upper state column density maps and error maps into lists
        log_N_upper_gs.append(log_N_upper_g)
        log_N_upper_g_errs.append(log_N_upper_g_err)

        # Plot moment 0 and upper state column density maps for each rung next to each other
    #     fig = plt.figure(figsize = (20, 10))
    #     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
    #     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
    #     plt.show()

        print(np.nanmin(log_N_upper_g))
        print(np.nanmax(log_N_upper_g))

        i += 1

    log_N_upper_gs = np.array(log_N_upper_gs)
    log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    return log_N_upper_gs, log_N_upper_g_errs
log_N_upper_gs, log_N_upper_g_errs = generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g)
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
## DO WE STILL NEED TO DO THIS?
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
# plt.show()
### NEW STUFF 11/10
temperature_final = (temperature*u.erg/constants.k_B).decompose()
log_intercept_final = np.log10(np.exp(ln_intercept)*ch3cn_partfunc(temperature_final))
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_intercept_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # used to be log_N_tot
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
# plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Export maps
# header = hdu[0].header
# fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
# fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
max_temp_coord = np.unravel_index(np.nanargmax(temperature_final), 
                                  temperature_final.shape)
print(max_temp_coord)
temp_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
intercept_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
print(temp_extract)
# ch3cn_E_U might be the wrong shape (might be 2D instead of 1D)
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['figure.background'] = 'white'
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['savefig.facecolor'] = 'white'
########################################################
# Started Logging At: 2021-12-01 14:41:40
########################################################
########################################################
# # Started Logging At: 2021-12-01 14:41:43
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['savefig.facecolor'] = 'white'
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
def retrieve_cube(results, freq_spw):
    '''Get methyl cyanide (target molecule) cube'''
    fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
    ch3cncube = SpectralCube.read(fn, format='casa_image')
    ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN',
                                                                                        fmin=fmin, 
                                                                                        fmax=fmax, 
                                                                                        catalog='JPL')
    # We're readying the partition function for use with temperature map later!
    ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
    ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data
    imshow(noise_map, origin='lower')
    return noise_map
noise_map = generate_noise_map()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        noise_map_int = noise_map*u.K*(channel_width)

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

        # Append upper state column density maps and error maps into lists
        log_N_upper_gs.append(log_N_upper_g)
        log_N_upper_g_errs.append(log_N_upper_g_err)

        # Plot moment 0 and upper state column density maps for each rung next to each other
    #     fig = plt.figure(figsize = (20, 10))
    #     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
    #     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
    #     plt.show()

        print(np.nanmin(log_N_upper_g))
        print(np.nanmax(log_N_upper_g))

        i += 1

    log_N_upper_gs = np.array(log_N_upper_gs)
    log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    return log_N_upper_gs, log_N_upper_g_errs
log_N_upper_gs, log_N_upper_g_errs = generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g)
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
## DO WE STILL NEED TO DO THIS?
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
# plt.show()
### NEW STUFF 11/10
temperature_final = (temperature*u.erg/constants.k_B).decompose()
log_intercept_final = np.log10(np.exp(ln_intercept)*ch3cn_partfunc(temperature_final))
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_intercept_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # used to be log_N_tot
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
# plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Export maps
# header = hdu[0].header
# fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
# fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
max_temp_coord = np.unravel_index(np.nanargmax(temperature_final), 
                                  temperature_final.shape)
print(max_temp_coord)
temp_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
intercept_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
print(temp_extract)
# ch3cn_E_U might be the wrong shape (might be 2D instead of 1D)
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['axes.facecolor'] = 'white'
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data
    imshow(noise_map, origin='lower')
    return noise_map
noise_map = generate_noise_map()
get_ipython().run_line_magic('matplotlib', 'inline')

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b3fd6e7ce10>
get_ipython().run_line_magic('matplotlib', 'inline')

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b3fd6f31828>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
header = hdu[0].header
fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
temperature[~np.isinf(temperature)].shape
# np.histogram(temperature[~np.isinf(temperature)])
n, bins, patches = plt.hist(temperature[~np.isinf(temperature)], bins=100)
plt.yscale('log')
# plt.xlim(-1000, 1000)
plt.xlabel("Fitted temperature value (K)")
plt.ylabel("Counts per bin")
plt.show()

print("Pixels fitted:",temperature.shape[0]*temperature.shape[1]) # How many temperature values did we calculate? (2D)
print("Above zero:",np.sum(temperature > 0)) # How many of the values are above zero?
print("Equal to zero:",np.sum(temperature == 0.0)) # How many values are equal to zero?
# test = [0, 1, 2, 3]
# test2 = test.copy()

all_comp_coord = np.unravel_index(np.argmax(ln_N_upper_gs_mask_sum), 
                                  temperature.shape)
ln_N_upper_gs_mask_sum.max()
#[Out]# 8
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParmas['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b3fe2150240>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
header = hdu[0].header
fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
temperature[~np.isinf(temperature)].shape
# np.histogram(temperature[~np.isinf(temperature)])
n, bins, patches = plt.hist(temperature[~np.isinf(temperature)], bins=100)
plt.yscale('log')
# plt.xlim(-1000, 1000)
plt.xlabel("Fitted temperature value (K)")
plt.ylabel("Counts per bin")
plt.show()

print("Pixels fitted:",temperature.shape[0]*temperature.shape[1]) # How many temperature values did we calculate? (2D)
print("Above zero:",np.sum(temperature > 0)) # How many of the values are above zero?
print("Equal to zero:",np.sum(temperature == 0.0)) # How many values are equal to zero?
ln_N_upper_gs_mask_sum.max()
#[Out]# 8
# test = [0, 1, 2, 3]
# test2 = test.copy()

all_comp_coord = np.unravel_index(np.argmax(ln_N_upper_gs_mask_sum), 
                                  temperature.shape)
# Plot rotational diagram for a pixel
max_temp_coord = np.unravel_index(np.nanargmax(temperature), 
                                  temperature.shape)
print(max_temp_coord)

# temperature_new = temperature.copy()
# temperature_new[max_temp_coord[0], max_temp_coord[1]] == np.nan
# max_temp_coord_new = np.unravel_index(np.nanargmax(temperature_new), 
#                                       temperature_new.shape)
# print(max_temp_coord_new)

example_N_upper = ln_N_upper_gs[:, all_comp_coord[0], all_comp_coord[1]]
example_E_upper = ch3cn_E_U

print(example_E_upper)

plt.style.use('dark_background')

fig = plt.figure(figsize=(7, 4))
plt.scatter(example_E_upper.value, example_N_upper)
plt.ylabel("$\log (N_u / g)$")
plt.xlabel("$E_u / k$ (K)")
# plt.ylim(9.30225, 9.30275)
plt.text(2.7e9, 25, f"T = {np.nanargmax(temperature)} K")
plt.text(2.7e9, 24.6, f"Coordinate = {max_temp_coord}")
plt.show()
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
fmin = 110330.34*u.MHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 110383.6*u.MHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# used to be 'CH3CN' as first argument
mod = np.nansum(mods,axis=0)*fillingfactor
mod_sp = pyspeckit.Spectrum(data = mod*u.K, xarr = sp.xarr.copy())
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['figure.facecolor'] = 'white'
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
def retrieve_cube(results, freq_spw):
    '''Get methyl cyanide (target molecule) cube'''
    fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
    ch3cncube = SpectralCube.read(fn, format='casa_image')
    ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN',
                                                                                        fmin=fmin, 
                                                                                        fmax=fmax, 
                                                                                        catalog='JPL')
    # We're readying the partition function for use with temperature map later!
    ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
    ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data
    imshow(noise_map, origin='lower')
    return noise_map
noise_map = generate_noise_map()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        noise_map_int = noise_map*u.K*(channel_width)

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

        # Append upper state column density maps and error maps into lists
        log_N_upper_gs.append(log_N_upper_g)
        log_N_upper_g_errs.append(log_N_upper_g_err)

        # Plot moment 0 and upper state column density maps for each rung next to each other
    #     fig = plt.figure(figsize = (20, 10))
    #     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
    #     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
    # #     plt.colorbar(mappable = im)
    #     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
    #     plt.show()

        print(np.nanmin(log_N_upper_g))
        print(np.nanmax(log_N_upper_g))

        i += 1

    log_N_upper_gs = np.array(log_N_upper_gs)
    log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    return log_N_upper_gs, log_N_upper_g_errs
log_N_upper_gs, log_N_upper_g_errs = generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g)
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
## DO WE STILL NEED TO DO THIS?
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
# plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
# plt.show()
### NEW STUFF 11/10
temperature_final = (temperature*u.erg/constants.k_B).decompose()
log_intercept_final = np.log10(np.exp(ln_intercept)*ch3cn_partfunc(temperature_final))
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_intercept_final.data, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # used to be log_N_tot
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
# plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Export maps
# header = hdu[0].header
# fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
# fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
max_temp_coord = np.unravel_index(np.nanargmax(temperature_final), 
                                  temperature_final.shape)
print(max_temp_coord)
temp_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
intercept_extract = temperature_final[max_temp_coord[0]+1, max_temp_coord[1]-10]
print(temp_extract)
# ch3cn_E_U might be the wrong shape (might be 2D instead of 1D)
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs, '.')
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b3fea8b8da0>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
#plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
#plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
########################################################
# Started Logging At: 2021-12-01 14:54:10
########################################################
########################################################
# # Started Logging At: 2021-12-01 14:54:12
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b95e32b4198>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
header = hdu[0].header
fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
temperature[~np.isinf(temperature)].shape
# np.histogram(temperature[~np.isinf(temperature)])
n, bins, patches = plt.hist(temperature[~np.isinf(temperature)], bins=100)
plt.yscale('log')
# plt.xlim(-1000, 1000)
plt.xlabel("Fitted temperature value (K)")
plt.ylabel("Counts per bin")
plt.show()

print("Pixels fitted:",temperature.shape[0]*temperature.shape[1]) # How many temperature values did we calculate? (2D)
print("Above zero:",np.sum(temperature > 0)) # How many of the values are above zero?
print("Equal to zero:",np.sum(temperature == 0.0)) # How many values are equal to zero?
ln_N_upper_gs_mask_sum.max()
#[Out]# 8
# test = [0, 1, 2, 3]
# test2 = test.copy()

all_comp_coord = np.unravel_index(np.argmax(ln_N_upper_gs_mask_sum), 
                                  temperature.shape)
# Plot rotational diagram for a pixel
max_temp_coord = np.unravel_index(np.nanargmax(temperature), 
                                  temperature.shape)
print(max_temp_coord)

# temperature_new = temperature.copy()
# temperature_new[max_temp_coord[0], max_temp_coord[1]] == np.nan
# max_temp_coord_new = np.unravel_index(np.nanargmax(temperature_new), 
#                                       temperature_new.shape)
# print(max_temp_coord_new)

example_N_upper = ln_N_upper_gs[:, all_comp_coord[0], all_comp_coord[1]]
example_E_upper = ch3cn_E_U

print(example_E_upper)

plt.style.use('dark_background')

fig = plt.figure(figsize=(7, 4))
plt.scatter(example_E_upper.value, example_N_upper)
plt.ylabel("$\log (N_u / g)$")
plt.xlabel("$E_u / k$ (K)")
# plt.ylim(9.30225, 9.30275)
plt.text(2.7e9, 25, f"T = {np.nanargmax(temperature)} K")
plt.text(2.7e9, 24.6, f"Coordinate = {max_temp_coord}")
plt.show()
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
fmin = 110330.34*u.MHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 110383.6*u.MHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# used to be 'CH3CN' as first argument
mod = np.nansum(mods,axis=0)*fillingfactor
mod_sp = pyspeckit.Spectrum(data = mod*u.K, xarr = sp.xarr.copy())
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b960235e940>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
########################################################
# Started Logging At: 2021-12-01 14:56:41
########################################################
########################################################
# # Started Logging At: 2021-12-01 14:56:42
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b60b29ea198>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
header = hdu[0].header
fits.PrimaryHDU(data=temperature, header=header).writeto('temperature_map_test.fits', overwrite=True)
fits.PrimaryHDU(data=log_N_tot, header=header).writeto('log_N_tot_test.fits', overwrite=True)
temperature
#[Out]# <Quantity [[-inf, -inf, -inf, ..., -inf, -inf, -inf],
#[Out]#            [-inf, -inf, -inf, ..., -inf, -inf, -inf],
#[Out]#            [-inf, -inf, -inf, ..., -inf, -inf, -inf],
#[Out]#            ...,
#[Out]#            [-inf, -inf, -inf, ..., -inf, -inf, -inf],
#[Out]#            [-inf, -inf, -inf, ..., -inf, -inf, -inf],
#[Out]#            [-inf, -inf, -inf, ..., -inf, -inf, -inf]] K>
np.nanmax(temperature)
#[Out]# <Quantity 2.4268058e+24 K>
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab
pylab.rcParams['figure.facecolor'] = 'white'

from pylab import imshow
from astropy.io import fits
hdu = fits.open('methyl_cyanide/template_noise.fits')
noise_map = hdu[0].data
imshow(noise_map, origin='lower') 
#[Out]# <matplotlib.image.AxesImage at 0x2b60bd955940>
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

# Get methyl cyanide (target molecule) cube
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
# Check angular area of pixels:
wcs = ch3cncube.wcs
wcs.proj_plane_pixel_area()
#[Out]# <Quantity 3.08641975e-09 deg2>
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg


# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# array([6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#        2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, shifted_line_freqs, 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
#ln_N_upper_gs = np.log(10**(log_N_upper_gs))
#ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
#ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
#ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ch3cn_E_U
#[Out]# <Quantity [3.81587925e+09, 2.88846564e+09, 2.10337313e+09, 1.46078157e+09,
#[Out]#            9.60836271e+08, 6.03650888e+08, 3.89307430e+08, 3.17854816e+08] K / J>
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
pl.imshow(mom0.value)
plt.imshow(mom0.value)
#[Out]# <matplotlib.image.AxesImage at 0x2b60be35da58>
plt.imshow(mom0.value)
plt.colorbar()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2b60b20df278>
plt.imshow(mom0.value)
plt.colorbar()
mom0.unit
#[Out]# Unit("K km / s")
shifted_line_freqs
#[Out]# <Projection [[nan, nan, nan, ..., nan, nan, nan],
#[Out]#              [nan, nan, nan, ..., nan, nan, nan],
#[Out]#              [nan, nan, nan, ..., nan, nan, nan],
#[Out]#              ...,
#[Out]#              [nan, nan, nan, ..., nan, nan, nan],
#[Out]#              [nan, nan, nan, ..., nan, nan, nan],
#[Out]#              [nan, nan, nan, ..., nan, nan, nan]] GHz>
plt.imshow(shifted_line_freqs.value)
plt.colorbar()
#[Out]# <matplotlib.colorbar.Colorbar at 0x2b60b2193940>
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, chc3n_freqs[i], 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.data/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
get_ipython().run_line_magic('pinfo2', 'nupper_of_kkms')
N_upper
#[Out]# <Quantity [[nan, nan, nan, ..., nan, nan, nan],
#[Out]#            [nan, nan, nan, ..., nan, nan, nan],
#[Out]#            [nan, nan, nan, ..., nan, nan, nan],
#[Out]#            ...,
#[Out]#            [nan, nan, nan, ..., nan, nan, nan],
#[Out]#            [nan, nan, nan, ..., nan, nan, nan],
#[Out]#            [nan, nan, nan, ..., nan, nan, nan]] 1 / cm2>
N_upper.data
#[Out]# <memory at 0x2b6070e34dc8>
N_upper.value
#[Out]# array([[nan, nan, nan, ..., nan, nan, nan],
#[Out]#        [nan, nan, nan, ..., nan, nan, nan],
#[Out]#        [nan, nan, nan, ..., nan, nan, nan],
#[Out]#        ...,
#[Out]#        [nan, nan, nan, ..., nan, nan, nan],
#[Out]#        [nan, nan, nan, ..., nan, nan, nan],
#[Out]#        [nan, nan, nan, ..., nan, nan, nan]])
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], 10**ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
np.nan_to_num
#[Out]# <function numpy.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None)>
get_ipython().run_line_magic('pinfo', 'np.nan_to_num')
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
# set the errors to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs
#[Out]# array([[[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]],
#[Out]# 
#[Out]#        [[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]],
#[Out]# 
#[Out]#        [[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]],
#[Out]# 
#[Out]#        ...,
#[Out]# 
#[Out]#        [[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]],
#[Out]# 
#[Out]#        [[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]],
#[Out]# 
#[Out]#        [[1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         ...,
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10],
#[Out]#         [1.e+10, 1.e+10, 1.e+10, ..., 1.e+10, 1.e+10, 1.e+10]]])
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = 1/ln_N_upper_g_errs**2
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
fitshape
#[Out]# (8, 262144)
ln_N_upper_g_errs.shape
#[Out]# (8, 512, 512)
ln_N_upper_g_errs.shape
#[Out]# (8, 512, 512)
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
weights.shape
#[Out]# (8, 262144)
ch3cn_E_U_1d.shape
#[Out]# (8, 2)
ffw.shape
#[Out]# (8, 8, 262144)
ffr.shape
#[Out]# (8, 262144)
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# weights = weights.reshape(fitshape)
ffw = (ffr * weights)[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# this works *if and only if* the weights are dependent on the Z-axis, but not the X- or Y- axis
# so we go back to basics: uniform weights (yikes)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
ffr.shape
#[Out]# (8, 262144)
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
ln_N_upper_gs.shape
#[Out]# (8, 512, 512)
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs[:,250,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
plt.imshow(log_N_upper_gs[0])
#[Out]# <matplotlib.image.AxesImage at 0x2b60bf117240>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
(ch3cn_E_U*u.erg/constants.k_B).decompose()
#[Out]# <Quantity [2.76383009e+25, 2.09210715e+25, 1.52346696e+25, 1.05803979e+25,
#[Out]#            6.95930878e+24, 4.37222558e+24, 2.81974224e+24, 2.30221306e+24] K2 s2 / (kg m2)>
(ch3cn_E_U*u.erg).decompose()
#[Out]# <Quantity [381.5879253 , 288.84656395, 210.33731322, 146.07815721,
#[Out]#             96.08362709,  60.36508878,  38.93074305,  31.78548157] K>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U*u.erg).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
print(ch3cn_E_U)

# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = ch3cn_E_U*u.erg/constants.k_B # Original is in erg
print(ch3cn_E_U)

# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A # Original is log_10(A_ij)
ch3cn_E_U = (ch3cn_E_U*u.erg/constants.k_B).decompose() # Original is in erg
print(ch3cn_E_U)

# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
# INSTEAD OF WHAT'S AT THE BOTTOM OF THIS CELL, USE THIS (without Tex set):
from lte_modeling_tools import get_molecular_parameters
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# We're readying the partition function for use with temperature map later!

from astropy import constants
ch3cn_A = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij)
ch3cn_E_U = (ch3cn_E_U*u.erg/constants.k_B).decompose() # Original is in erg
print(ch3cn_E_U)

# # from astroquery.splatalogue import Splatalogue
# ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
#                                 show_upper_degeneracy=True, show_qn_code=True)
# ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
# ch3cntbl = ch3cntbl[::-1]
# # ch3cn_freqs = ch3cntbl['Freq-GHz(rest frame,redshifted)']
# # ch3cn_A = 10**ch3cntbl['Log<sub>10</sub> (A<sub>ij</sub>)']
# # ch3cn_g = ch3cntbl['Upper State Degeneracy']
# # ch3cn_E_U = ch3cntbl['E_U (K)'][:,None,None] # This is already divided by the Boltzmann constant, I think?
# ch3cntbl.show_in_notebook()
ch3cn_A
#[Out]# <Quantity [6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#            2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04] 1 / s>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U*u.erg).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U*).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
ln_N_upper_g_errs.shape
#[Out]# (8, 512, 512)
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# this works *if and only if* the weights are dependent on the Z-axis, but not the X- or Y- axis
# so we go back to basics: uniform weights (yikes)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
ch3cn_A
#[Out]# <Quantity [6.27603981e-05, 1.17246369e-04, 1.63403281e-04, 2.01236068e-04,
#[Out]#            2.30651121e-04, 2.51719255e-04, 2.64307995e-04, 2.68580345e-04] 1 / s>
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
import matplotlib.pyplot as plt

log_N_upper_gs = []
log_N_upper_g_errs = []

i = 0
for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]

    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()
    mom1 = masked_cube.moment1()
    # Propagate error on integrated intensity
    noise_map_int = noise_map*u.K*(channel_width)

    # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
    #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

    # Calculate upper state column density from integrated line intensity (moment 0 map)
    N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], ch3cn_A[i])
#     print(np.nanmean(N_upper))
    log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .data?
#     print(np.nanmean(log_N_upper_g))
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, shifted_line_freqs, ch3cn_A[i])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant?

    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)
    
    # Plot moment 0 and upper state column density maps for each rung next to each other
#     fig = plt.figure(figsize = (20, 10))
#     im = plt.subplot(1,2,1).imshow(mom0.value, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Moment 0 Map")
#     im = plt.subplot(1,2,2).imshow(log_N_upper_g, origin='lower', cmap='magma')
# #     plt.colorbar(mappable = im)
#     plt.title(f"Methyl Cyanide $j = 8$, $k = {i}$ Upper State Column Density Map")
#     plt.show()

    print(np.nanmin(log_N_upper_g)) # These values seem a little low
    print(np.nanmax(log_N_upper_g))
    
    i += 1

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
masked_cube.max()

masked_cube_test = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
np.nanmax(masked_cube_test.moment0())
#[Out]# <Projection 8.10930887 K km / s>
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
# set the errors to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
## Some stuff that didn't work:
# valid_map = np.isfinite(log_N_upper_gs)
# valid_map.shape
# # Build up E_U matrix by hand, to be masked (to go with second fitting attempt)
# E_U_matrix = []
# i = 0
# for i in range(0, 8):
#     new_map = np.ones((512, 512))
#     new_map *= ch3cn_E_U[i,0,0]
#     E_U_matrix.append(new_map)
#     i += 1
# E_U_matrix = np.array(E_U_matrix)
# np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
# log_N_upper_gs
# np.moveaxis(np.tile(1/(noise_map**2)[:,:,None], 8), 2, 0).shape
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
ch3cn_E_U
#[Out]# <Quantity [381.5879253 , 288.84656395, 210.33731322, 146.07815721,
#[Out]#             96.08362709,  60.36508878,  38.93074305,  31.78548157] K>
plt.imshow(log_N_upper_gs[0])
#[Out]# <matplotlib.image.AxesImage at 0x2b60bf357ac8>
(ch3cn_E_U*u.erg).decompose()
#[Out]# <Quantity [3.81587925e-05, 2.88846564e-05, 2.10337313e-05, 1.46078157e-05,
#[Out]#            9.60836271e-06, 6.03650888e-06, 3.89307430e-06, 3.17854816e-06] K kg m2 / s2>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
ln_N_upper_g_errs.shape
#[Out]# (8, 512, 512)
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U, np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# this works *if and only if* the weights are dependent on the Z-axis, but not the X- or Y- axis
# so we go back to basics: uniform weights (yikes)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
ffr.shape
#[Out]# (8, 262144)
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U[::-1]).decompose(), log_N_upper_gs[:,200,300], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
plt.imshow(log_N_upper_gs[5])
#[Out]# <matplotlib.image.AxesImage at 0x2b60beeb8588>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U[::-1]).decompose(), log_N_upper_gs[:,197,298], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U[::-1], np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# this works *if and only if* the weights are dependent on the Z-axis, but not the X- or Y- axis
# so we go back to basics: uniform weights (yikes)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
ffr.shape
#[Out]# (8, 262144)
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma', vmax=500); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma', vmax=500, vmin=0); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
# sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # this does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 10 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
ch3cn_E_U
#[Out]# <Quantity [381.5879253 , 288.84656395, 210.33731322, 146.07815721,
#[Out]#             96.08362709,  60.36508878,  38.93074305,  31.78548157] K>
plt.imshow(log_N_upper_gs[5])
#[Out]# <matplotlib.image.AxesImage at 0x2b60c8a5cb70>
(ch3cn_E_U*u.erg).decompose()
#[Out]# <Quantity [3.81587925e-05, 2.88846564e-05, 2.10337313e-05, 1.46078157e-05,
#[Out]#            9.60836271e-06, 6.03650888e-06, 3.89307430e-06, 3.17854816e-06] K kg m2 / s2>
import matplotlib.pyplot as plt
plt.plot((ch3cn_E_U[::-1]).decompose(), log_N_upper_gs[:,197,298], '.')
#plt.plot((ch3cn_E_U*u.erg/constants.k_B).decompose(), np.log10(np.exp(fitted_m*(ch3cn_E_U)+fitted_b)))
plt.show()
ln_N_upper_g_errs.shape
#[Out]# (8, 512, 512)
# Original fitting attempt
# ch3cn_E_U_1d = np.array([ch3cn_E_U[:,0,0], np.ones(ch3cn_E_U.shape[0])]).T
ch3cn_E_U_1d = np.array([ch3cn_E_U[::-1], np.ones(ch3cn_E_U.shape[0])]).T
fitshape = ln_N_upper_gs.shape[0], np.product(ln_N_upper_gs.shape[1:])
ffr = ln_N_upper_gs.reshape(fitshape)
# ffw = ffr * 1 # Ignoring weights for now 
# zzw = ch3cn_E_U_1d * 1 # Ignoring weights for now

# Do weights stuff: 
#weights = np.mean(ln_N_upper_g_errs, axis=(1,2))
weights = (1/ln_N_upper_g_errs**2).reshape(fitshape)
# this works *if and only if* the weights are dependent on the Z-axis, but not the X- or Y- axis
# so we go back to basics: uniform weights (yikes)
weights = np.ones(8)
# weights = weights.reshape(fitshape)
ffw = ffr * weights[:,None]
zzw = ch3cn_E_U_1d * weights[:,None]

print(f"ffw.shape={ffw.shape}, zz1d.shape={ch3cn_E_U_1d.shape}")
fitted, residuals, rank, singular = np.linalg.lstsq(zzw, ffw, rcond=None)
print(f"rank={rank}, singular={singular}")

# Extract fit values
fitted_m = fitted[0].reshape(ln_N_upper_gs.shape[1:]) # not sure what fitted[0] is vs. fitted[1]
fitted_b = fitted[1].reshape(ln_N_upper_gs.shape[1:]) # m or b? we just don't know
fitted_resid = residuals.reshape(ln_N_upper_gs.shape[1:])

temperature = -1./fitted_m # okay, we are assuming that this is the correct interpretation of m vs. b
temperature = (temperature*u.erg/constants.k_B).decompose()
ln_intercept = fitted_b
# log_intercept = np.log10(np.exp(ln_intercept)) # Log base 10 of (e to the [natural log of N])

## This is wrong:
# temperature_alt = -1./fitted_b
# ln_N_alt = fitted_m
# log_N_alt = np.log10(np.exp(ln_N_alt))
ffr.shape
#[Out]# (8, 262144)
# log_N[196,291]
# N_tot = intercept * Z, so ln_N_tot = ln_intercept + ln_Z
Z = ch3cn_partfunc(temperature)
ln_Z = np.log(Z)
ln_N_tot = ln_intercept + ln_Z
log_N_tot = np.log10(np.exp(ln_N_tot)) # Log base 10 of (e to the [natural log of N])
# Plot results
fig = plt.figure(figsize = (15, 15))
im = plt.subplot(2,2,1).imshow(temperature.value, origin='lower', cmap='magma', vmax=500, vmin=0); plt.colorbar(mappable = im, fraction=0.046, pad=0.04) # , vmin=0, vmax=100
plt.title("Temperature (K)")
im = plt.subplot(2,2,2).imshow(log_N_tot, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
plt.title("Column density ($\log_{10}(N_{tot})$)")
# im = plt.subplot(2,2,3).imshow(temperature_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Temperature (K, alternate)")
# im = plt.subplot(2,2,4).imshow(log_N_alt, origin='lower', cmap='magma'); plt.colorbar(mappable = im, fraction=0.046, pad=0.04)
# plt.title("Column density (units?, alternate)")
# im = plt.subplot(2,2,4).imshow(fitted_resid, origin='lower', cmap='magma'); plt.colorbar(mappable = im)
# plt.title("Fitted residuals")
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);
plt.savefig(f'/blue/adamginsburg/abulatek/brick/first_results/prelim_temp_map.pdf')
plt.show()

# Making colorbars smaller
# https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
########################################################
# Started Logging At: 2021-12-01 15:44:18
########################################################
########################################################
# # Started Logging At: 2021-12-01 15:44:21
########################################################
########################################################
# Started Logging At: 2021-12-01 15:44:27
########################################################
########################################################
# # Started Logging At: 2021-12-01 15:44:27
########################################################
