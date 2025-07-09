########################################################
# Started Logging At: 2022-02-16 16:18:08
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:18:10
########################################################
########################################################
# Started Logging At: 2022-02-16 16:19:38
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:19:39
########################################################
########################################################
# Started Logging At: 2022-02-16 16:19:41
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:19:41
########################################################
########################################################
# Started Logging At: 2022-02-16 16:19:58
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:19:58
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')

import pylab as pl
# pl.style.use('dark_background')
# pl.rcParams['figure.facecolor'] = 'w'

# plt.rc('text', usetex = True) # Use LaTeX font in plots
# plt.rc('font', family = 'serif')
# plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
# pl.rcParams.update({'font.size': 12})
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 200
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
    ch3cn_A = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij)
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    #ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution
                                # is not implemented yet
    imshow(noise_map, origin='lower')
    return noise_map
noise_map = generate_noise_map()
# ch3cncube.mad_std(axis = 0)
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
# masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
# plt.imshow(masked_cube[0])
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
#         primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb', format='casa_image')
#         masked_cube = masked_cube/primary_beam # Correct for effect of primary beam
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)

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
# Set the errors for NaN values to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
# Don't fit missing data (SHOULD INCORP UPPER LIMITS)
pix_x, pix_y = 319, 220
col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pix_y,pix_x]) # ignores NaNs
col_density_zero_mask = ln_N_upper_gs[:,pix_y,pix_x] != 0
fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
fit_col_densities = ln_N_upper_gs[:,pix_y,pix_x][col_density_zero_mask]
fit_col_densities_errs = ln_N_upper_g_errs[:,pix_y,pix_x][col_density_zero_mask]

print(fit_energies)
print(fit_col_densities)
print(fit_col_densities_errs)
# do the simple linear fit
from scipy.optimize import curve_fit

def linear(x, m, b):
    return m*x + b

# CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()

# guess = np.array([1e-15, 100])
# Fit without errorbars:
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# Fit with errorbars: 
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, 
                       sigma = fit_col_densities_errs) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
fig = plt.figure(dpi = 200) # I shouldn't have to do this every time
plt.errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)),
             yerr = np.log10(fit_col_densities_errs), fmt = '.') # should these both be log10'ed?
plt.plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)))
plt.show()
from astropy.visualization import quantity_support
quantity_support()
#[Out]# <astropy.visualization.units.quantity_support.<locals>.MplQuantityConverter at 0x2acbe7203040>
noise_map = generate_noise_map()
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution
                                # is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
from astropy.visualization import quantity_support
quantity_support()
#[Out]# <astropy.visualization.units.quantity_support.<locals>.MplQuantityConverter at 0x2acbe721ea90>
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution
                                # is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
from astropy.visualization import quantity_support
quantity_support()
#[Out]# <astropy.visualization.units.quantity_support.<locals>.MplQuantityConverter at 0x2acbe723fb20>
noise_map = generate_noise_map()
# ch3cncube.mad_std(axis = 0)
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
# masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
# plt.imshow(masked_cube[0])
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
#         primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb', format='casa_image')
#         masked_cube = masked_cube/primary_beam # Correct for effect of primary beam
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)

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
# Set the errors for NaN values to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
# Don't fit missing data (SHOULD INCORP UPPER LIMITS)
pix_x, pix_y = 319, 220
col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pix_y,pix_x]) # ignores NaNs
col_density_zero_mask = ln_N_upper_gs[:,pix_y,pix_x] != 0
fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
fit_col_densities = ln_N_upper_gs[:,pix_y,pix_x][col_density_zero_mask]
fit_col_densities_errs = ln_N_upper_g_errs[:,pix_y,pix_x][col_density_zero_mask]

print(fit_energies)
print(fit_col_densities)
print(fit_col_densities_errs)
# do the simple linear fit
from scipy.optimize import curve_fit

def linear(x, m, b):
    return m*x + b

# CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()

# guess = np.array([1e-15, 100])
# Fit without errorbars:
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# Fit with errorbars: 
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, 
                       sigma = fit_col_densities_errs) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
fig = plt.figure(dpi = 200) # I shouldn't have to do this every time
plt.errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)),
             yerr = np.log10(fit_col_densities_errs), fmt = '.') # should these both be log10'ed?
plt.plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)))
plt.show()
# fit_energies_converted[3]
# fit_col_densities[3]
# channel_width
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
channel_width = np.diff(masked_cube.spectral_axis)[0]

print("Noise map:",noise_map[pix_y,pix_x])

number_of_pixels = masked_cube.mask.include().sum(axis=0)[pix_y,pix_x]

noise_map_int = noise_map[pix_y,pix_x]*channel_width*np.sqrt(number_of_pixels)
print("Noise map integrated:",noise_map_int)

mom0 = masked_cube.moment0()
print("mom0 value at pixel:",mom0[pix_y,pix_x])

N_upper = nupper_of_kkms(mom0[pix_y,pix_x], ch3cn_freqs[3], ch3cn_A[3])
print("N_upper:",N_upper)
N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[3], ch3cn_A[3])
print("Error on N_upper:",N_upper_err)
log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.))
print("Error on log of N_upper:",log_N_upper_g_err)
# # extract the values
# temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
# # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
# total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

# print(f"Temp: {temp}")
# print(f"log10(Total column density): {total_col_density}")
from scipy.optimize import curve_fit
import tqdm
print(tqdm.__version__)
from tqdm.notebook import tqdm
########################################################
# Started Logging At: 2022-02-16 16:24:07
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:24:08
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')

import pylab as pl
# pl.style.use('dark_background')
# pl.rcParams['figure.facecolor'] = 'w'

# plt.rc('text', usetex = True) # Use LaTeX font in plots
# plt.rc('font', family = 'serif')
# plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
# pl.rcParams.update({'font.size': 12})
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 200
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
    ch3cn_A = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij)
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    #ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution
                                # is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
from astropy.visualization import quantity_support
quantity_support()
#[Out]# <astropy.visualization.units.quantity_support.<locals>.MplQuantityConverter at 0x2b21ff815940>
noise_map = generate_noise_map()
# ch3cncube.mad_std(axis = 0)
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
# masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
# plt.imshow(masked_cube[0])
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
#         primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb', format='casa_image')
#         masked_cube = masked_cube/primary_beam # Correct for effect of primary beam
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)

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
# Set the errors for NaN values to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
# Don't fit missing data (SHOULD INCORP UPPER LIMITS)
pix_x, pix_y = 319, 220
col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pix_y,pix_x]) # ignores NaNs
col_density_zero_mask = ln_N_upper_gs[:,pix_y,pix_x] != 0
fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
fit_col_densities = ln_N_upper_gs[:,pix_y,pix_x][col_density_zero_mask]
fit_col_densities_errs = ln_N_upper_g_errs[:,pix_y,pix_x][col_density_zero_mask]

print(fit_energies)
print(fit_col_densities)
print(fit_col_densities_errs)
# do the simple linear fit
from scipy.optimize import curve_fit

def linear(x, m, b):
    return m*x + b

# CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()

# guess = np.array([1e-15, 100])
# Fit without errorbars:
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# Fit with errorbars: 
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, 
                       sigma = fit_col_densities_errs) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
fig = plt.figure(dpi = 200) # I shouldn't have to do this every time
plt.errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)),
             yerr = np.log10(fit_col_densities_errs), fmt = '.') # should these both be log10'ed?
plt.plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)))
plt.show()
# fit_energies_converted[3]
# fit_col_densities[3]
# channel_width
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
channel_width = np.diff(masked_cube.spectral_axis)[0]

print("Noise map:",noise_map[pix_y,pix_x])

number_of_pixels = masked_cube.mask.include().sum(axis=0)[pix_y,pix_x]

noise_map_int = noise_map[pix_y,pix_x]*channel_width*np.sqrt(number_of_pixels)
print("Noise map integrated:",noise_map_int)

mom0 = masked_cube.moment0()
print("mom0 value at pixel:",mom0[pix_y,pix_x])

N_upper = nupper_of_kkms(mom0[pix_y,pix_x], ch3cn_freqs[3], ch3cn_A[3])
print("N_upper:",N_upper)
N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[3], ch3cn_A[3])
print("Error on N_upper:",N_upper_err)
log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.))
print("Error on log of N_upper:",log_N_upper_g_err)
# # extract the values
# temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
# # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
# total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

# print(f"Temp: {temp}")
# print(f"log10(Total column density): {total_col_density}")
from scipy.optimize import curve_fit
import tqdm
print(tqdm.__version__)
from tqdm.notebook import tqdm
def rotational_diagram(pixel_x, pixel_y, plot = False, save = False):
    #  Don't fit missing data (SHOULD INCORP UPPER LIMITS)
    col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pixel_y,pixel_x]) # ignores NaNs
    col_density_zero_mask = ln_N_upper_gs[:,pixel_y,pixel_x] != 0
    fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
    fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x][col_density_zero_mask]
    fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x][col_density_zero_mask]
#     print(fit_energies, fit_col_densities)
    
    # do the simple linear fit
    # CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
    fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()
    # guess = np.array([1e-15, 100])
    popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, sigma = fit_col_densities_errs)
#     print(popt)
    slope, intercept = popt[0], popt[1]
    
    # extract the values
    temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
    # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
    total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

    if plot == True:
        # x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
        plt.errorbar(fit_energies_converted, fit_col_densities, yerr = fit_col_densities_errs, 
                     fmt = 'o', color = 'tab:blue') 
        plt.plot(fit_energies_converted, slope*fit_energies_converted.value+intercept, 
                 label = f"$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {total_col_density:.2f}",
                 color = 'tab:orange') # Why were we exp'ing the y axis???
    #     plt.yscale('log')
        plt.title(f"Rotational diagram for ({pixel_x}, {pixel_y})")
        plt.xlabel(f"Upper state energy [??? K but how get units?]")
        plt.ylabel(f"Column density [??? cm^{-2} but how get units?]")
        plt.legend()
        if save == True:
            plt.savefig("../figures/rot_diagram.pdf", dpi = 200, facecolor='w', edgecolor='w', bbox_inches='tight')
        plt.show()

        # print extracted values
        print(f"Temp: {temp:.5f}")
        print(f"log10(Total column density): {total_col_density:.5f}")

    return temp.value, total_col_density

def linear(x, m, b):
    return m*x + b
energy_K = (ch3cn_E_U[::-1]*u.erg/constants.k_B).decompose()
temp_map = np.zeros(ln_N_upper_gs.shape[1:])
total_col_density_map = np.zeros(ln_N_upper_gs.shape[1:])
for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
    col_density = ln_N_upper_gs[:,index[0],index[1]]
    col_density_errs = ln_N_upper_g_errs[:,index[0],index[1]]
    # Create column density mask, accounting for NaNs and zeros
    col_density_mask = (col_density != 0) & np.isfinite(col_density)
    if not col_density_mask.sum() > 1: # Are there unmasked values? If no, skip index
        continue
    temp_map[index], total_col_density_map[index] = rotational_diagram(index[1], index[0]) # Should be x, y order for function
# # extract the values
# temp_map = (-1./slope_map)*u.K
# # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
# total_col_density_map = np.log10(np.exp(intercept_map)*ch3cn_partfunc(temp_map))
# Set up colormap
# rock = sns.color_palette("rocket", as_cmap=True)
fig = plt.figure(dpi = 200)
cm = pl.matplotlib.cm.plasma.copy()
cm.set_under('w') # Make sure the "zero" color is white

# Get WCS coordinates from an (arbitrary?) masked cube
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
mom0 = masked_cube.moment0()

plt.subplot(111,  projection = mom0.wcs)
plt.imshow(temp_map, vmax=200, vmin=0.001, cmap = cm, origin='lower')
plt.tick_params(direction = 'in')
plt.title('Temperature (K)')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar()
plt.savefig("../figures/temp_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
# mako_r = sns.color_palette("mako_r", as_cmap=True)

fig = plt.figure(dpi = 200)

# Set up normalization
from astropy.visualization import imshow_norm, ManualInterval, SqrtStretch, SinhStretch

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w') # Make sure the "zero" color is white
# pl.imshow(total_col_density_map, cmap = cm, norm = norm, origin='lower')

plt.subplot(111,  projection = mom0.wcs)
im, norm = imshow_norm(total_col_density_map, origin='lower', cmap = cm, 
                       interval=ManualInterval(vmin = 10, vmax = np.max(total_col_density_map)), 
                       stretch=SinhStretch())
plt.tick_params(direction = 'in')
plt.title('Total Column Density ($\log_{10}$)')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar(mappable=im)
plt.savefig("../figures/coldens_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
fig = plt.figure(dpi = 200) # I DEFINITELY shouldn't have to run this again
rotational_diagram(319, 220, plot = True, save = True)
#[Out]# (38.324439210446414, 13.716505349556444)
pixel_x, pixel_y = 319, 220

### Get mom0 for each component at this pixel (and noise)
log_N_upper_gs = []
log_N_upper_g_errs = []
for comp in np.arange(0, 8, 1):
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{comp}_masked.fits', format='fits')
    channel_width = np.diff(masked_cube.spectral_axis)[0]
    # Calculate moment 0 and moment 1 maps of cube
    mom0 = masked_cube.moment0()[pixel_y, pixel_x]
    # Propagate error on integrated intensity
    number_of_pixels = masked_cube.mask.include().sum(axis=0)
    noise_map_int = noise_map[pixel_y, pixel_x]*channel_width*np.sqrt(number_of_pixels)
        # Why were we not multiplying by the channel width???
    ### Convert to N_upper (propagate error)
    N_upper = nupper_of_kkms(mom0, ch3cn_freqs[comp], ch3cn_A[comp])
    log_N_upper_g = np.log10(N_upper.value/ch3cn_g[comp]) # Shouldn't have to do .value?
    # Propagate error on upper state column density
    N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[comp], ch3cn_A[comp])
    log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)
    # Append upper state column density maps and error maps into lists
    log_N_upper_gs.append(log_N_upper_g)
    log_N_upper_g_errs.append(log_N_upper_g_err)

log_N_upper_gs = np.array(log_N_upper_gs)
log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    
### Filter out low SNR data
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))
# Replace all NaNs with 0s
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
# Set the errors for NaN values to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)
ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities
ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these

### Perform fit
col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pixel_y,pixel_x]) # ignores NaNs
col_density_zero_mask = ln_N_upper_gs[:,pixel_y,pixel_x] != 0
fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x][col_density_zero_mask]
fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x][col_density_zero_mask]
# CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, sigma = fit_col_densities_errs)
slope, intercept = popt[0], popt[1]
temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

### Plot rotational diagram with fit, print out fit temp and N_tot
plt.errorbar(fit_energies_converted, fit_col_densities, yerr = fit_col_densities_errs, 
             fmt = 'o', color = 'tab:blue') 
plt.plot(fit_energies_converted, slope*fit_energies_converted.value+intercept, 
         label = f"$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {total_col_density:.2f}",
         color = 'tab:orange') # Why were we exp'ing the y axis???
# plt.yscale('log')
plt.title(f"Rotational diagram for ({pixel_x}, {pixel_y})")
plt.xlabel(f"Upper state energy [??? K but how get units?]")
plt.ylabel(f"Column density [??? cm^{-2} but how get units?]")
plt.legend()
plt.show()

### Show data spectrum on synthetic spectrum with that temp and N_tot
fmin = 147.1*u.GHz
fmax = 147.2*u.GHz
sp_axis = np.linspace(fmin, fmax, 1000)
fillingfactor = 1
offset = 0
species = 'CH3CN'
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
# temp = temp
N_tot = total_column_density
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")
plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
# rotational_diagram(250, 300)
# rotational_diagram(255, 300)
# rotational_diagram(343, 278)
# rotational_diagram(345, 280)
# rotational_diagram(325, 281)
## Step 8: Generate synthetic spectrum for a pixel, and compare to original spectrum

# # Retrieve and plot original spectrum for one of these pixels (takes a while to run)
# test_pixel_x, test_pixel_y = 250, 300
# spectrum = ch3cncube[:,test_pixel_y,test_pixel_x]

# # Plot original spectrum for one of these pixels
# spectrum.quicklook()

# # Print temp/col density inputs from fit for this pixel
# rotational_diagram(test_pixel_x, test_pixel_y)

# %%capture
# # Generate synthetic spectrum for this cube, making sure to use temp/col density inputs from fit
# %run /blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py

# pl.plot(sp.xarr, sp.data)

# pl.plot(spectrum.value-np.array(sp.data))
# pl.title("Residuals (data minus synthetic), 1st attempt")
# pl.xlabel("Spectral units :)")
# pl.ylabel('Data minus model (maybe not same units?)')
from astropy import units as u

fmin = 147.1*u.GHz
fmax = 147.2*u.GHz
import numpy as np

sp_axis = np.linspace(fmin, fmax, 1000)
fillingfactor = 1
offset = 0
species = 'CH3CN'
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule

freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
temp = 38.32444*u.K
N_tot = (10**(13.71651))*u.cm**-2
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
pixel_x, pixel_y = 319, 220
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")
plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
fig = plt.figure(dpi = 200) # I DEFINITELY shouldn't have to run this again
temperature, total_col_density = rotational_diagram(319, 220, plot = True, save = True)
########################################################
# Started Logging At: 2022-02-16 16:28:00
########################################################
########################################################
# # Started Logging At: 2022-02-16 16:28:00
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')

import pylab as pl
# pl.style.use('dark_background')
# pl.rcParams['figure.facecolor'] = 'w'

# plt.rc('text', usetex = True) # Use LaTeX font in plots
# plt.rc('font', family = 'serif')
# plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
# pl.rcParams.update({'font.size': 12})
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 200
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
    ch3cn_A = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij)
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    #ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = retrieve_cube(results, freq_spw)
from pylab import imshow
from astropy.io import fits
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution
                                # is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
from astropy.visualization import quantity_support
quantity_support()
#[Out]# <astropy.visualization.units.quantity_support.<locals>.MplQuantityConverter at 0x2b3d02fa9310>
noise_map = generate_noise_map()
# ch3cncube.mad_std(axis = 0)
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
# masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
# plt.imshow(masked_cube[0])
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []

    i = 0
    for i in range(0, 8):

        # Import masked cube and get channel width
        masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
#         primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb', format='casa_image')
#         masked_cube = masked_cube/primary_beam # Correct for effect of primary beam
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Calculate moment 0 and moment 1 maps of cube
        mom0 = masked_cube.moment0()
        mom1 = masked_cube.moment1()
        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[i]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[i], ch3cn_A[i])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)

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
# Set the errors for NaN values to be huge so these are ignored in the fit
ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # now officially applying to data! *****
ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # errors big = try to ignore these
# Don't fit missing data (SHOULD INCORP UPPER LIMITS)
pix_x, pix_y = 319, 220
col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pix_y,pix_x]) # ignores NaNs
col_density_zero_mask = ln_N_upper_gs[:,pix_y,pix_x] != 0
fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
fit_col_densities = ln_N_upper_gs[:,pix_y,pix_x][col_density_zero_mask]
fit_col_densities_errs = ln_N_upper_g_errs[:,pix_y,pix_x][col_density_zero_mask]

print(fit_energies)
print(fit_col_densities)
print(fit_col_densities_errs)
# do the simple linear fit
from scipy.optimize import curve_fit

def linear(x, m, b):
    return m*x + b

# CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()

# guess = np.array([1e-15, 100])
# Fit without errorbars:
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# Fit with errorbars: 
popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, 
                       sigma = fit_col_densities_errs) #,  guess
print(popt)
slope, intercept = popt[0], popt[1]
# x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
fig = plt.figure(dpi = 200) # I shouldn't have to do this every time
plt.errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)),
             yerr = np.log10(fit_col_densities_errs), fmt = '.') # should these both be log10'ed?
plt.plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)))
plt.show()
# fit_energies_converted[3]
# fit_col_densities[3]
# channel_width
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_3_masked.fits', format='fits')
channel_width = np.diff(masked_cube.spectral_axis)[0]

print("Noise map:",noise_map[pix_y,pix_x])

number_of_pixels = masked_cube.mask.include().sum(axis=0)[pix_y,pix_x]

noise_map_int = noise_map[pix_y,pix_x]*channel_width*np.sqrt(number_of_pixels)
print("Noise map integrated:",noise_map_int)

mom0 = masked_cube.moment0()
print("mom0 value at pixel:",mom0[pix_y,pix_x])

N_upper = nupper_of_kkms(mom0[pix_y,pix_x], ch3cn_freqs[3], ch3cn_A[3])
print("N_upper:",N_upper)
N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[3], ch3cn_A[3])
print("Error on N_upper:",N_upper_err)
log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.))
print("Error on log of N_upper:",log_N_upper_g_err)
# # extract the values
# temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
# # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
# total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

# print(f"Temp: {temp}")
# print(f"log10(Total column density): {total_col_density}")
from scipy.optimize import curve_fit
import tqdm
print(tqdm.__version__)
from tqdm.notebook import tqdm
def rotational_diagram(pixel_x, pixel_y, plot = False, save = False):
    #  Don't fit missing data (SHOULD INCORP UPPER LIMITS)
    col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pixel_y,pixel_x]) # ignores NaNs
    col_density_zero_mask = ln_N_upper_gs[:,pixel_y,pixel_x] != 0
    fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
    fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x][col_density_zero_mask]
    fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x][col_density_zero_mask]
#     print(fit_energies, fit_col_densities)
    
    # do the simple linear fit
    # CONVERT TO K FROM ERGS TO DO THE FIT BC ERGS HAD SMALL NUMBERS
    fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()
    # guess = np.array([1e-15, 100])
    popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, sigma = fit_col_densities_errs)
#     print(popt)
    slope, intercept = popt[0], popt[1]
    
    # extract the values
    temp = (-1./slope)*u.K # okay, we are assuming that this is the correct interpretation of m vs. b
    # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
    total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

    if plot == True:
        # x_axis = (fit_energies*u.erg/constants.k_B).decompose() # NOW X AXIS IS IN KELVIN
        plt.errorbar(fit_energies_converted, fit_col_densities, yerr = fit_col_densities_errs, 
                     fmt = 'o', color = 'tab:blue') 
        plt.plot(fit_energies_converted, slope*fit_energies_converted.value+intercept, 
                 label = f"$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {total_col_density:.2f}",
                 color = 'tab:orange') # Why were we exp'ing the y axis???
    #     plt.yscale('log')
        plt.title(f"Rotational diagram for ({pixel_x}, {pixel_y})")
        plt.xlabel(f"Upper state energy [??? K but how get units?]")
        plt.ylabel(f"Column density [??? cm^{-2} but how get units?]")
        plt.legend()
        if save == True:
            plt.savefig("../figures/rot_diagram.pdf", dpi = 200, facecolor='w', edgecolor='w', bbox_inches='tight')
        plt.show()

        # print extracted values
        print(f"Temp: {temp:.5f}")
        print(f"log10(Total column density): {total_col_density:.5f}")

    return temp.value, total_col_density

def linear(x, m, b):
    return m*x + b
energy_K = (ch3cn_E_U[::-1]*u.erg/constants.k_B).decompose()
temp_map = np.zeros(ln_N_upper_gs.shape[1:])
total_col_density_map = np.zeros(ln_N_upper_gs.shape[1:])
for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
    col_density = ln_N_upper_gs[:,index[0],index[1]]
    col_density_errs = ln_N_upper_g_errs[:,index[0],index[1]]
    # Create column density mask, accounting for NaNs and zeros
    col_density_mask = (col_density != 0) & np.isfinite(col_density)
    if not col_density_mask.sum() > 1: # Are there unmasked values? If no, skip index
        continue
    temp_map[index], total_col_density_map[index] = rotational_diagram(index[1], index[0]) # Should be x, y order for function
# # extract the values
# temp_map = (-1./slope_map)*u.K
# # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BC FIT ALREADY DIVIDED BY K_B
# total_col_density_map = np.log10(np.exp(intercept_map)*ch3cn_partfunc(temp_map))
# Set up colormap
# rock = sns.color_palette("rocket", as_cmap=True)
fig = plt.figure(dpi = 200)
cm = pl.matplotlib.cm.plasma.copy()
cm.set_under('w') # Make sure the "zero" color is white

# Get WCS coordinates from an (arbitrary?) masked cube
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
mom0 = masked_cube.moment0()

plt.subplot(111,  projection = mom0.wcs)
plt.imshow(temp_map, vmax=200, vmin=0.001, cmap = cm, origin='lower')
plt.tick_params(direction = 'in')
plt.title('Temperature (K)')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar()
plt.savefig("../figures/temp_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
# mako_r = sns.color_palette("mako_r", as_cmap=True)

fig = plt.figure(dpi = 200)

# Set up normalization
from astropy.visualization import imshow_norm, ManualInterval, SqrtStretch, SinhStretch

cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w') # Make sure the "zero" color is white
# pl.imshow(total_col_density_map, cmap = cm, norm = norm, origin='lower')

plt.subplot(111,  projection = mom0.wcs)
im, norm = imshow_norm(total_col_density_map, origin='lower', cmap = cm, 
                       interval=ManualInterval(vmin = 10, vmax = np.max(total_col_density_map)), 
                       stretch=SinhStretch())
plt.tick_params(direction = 'in')
plt.title('Total Column Density ($\log_{10}$)')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar(mappable=im)
plt.savefig("../figures/coldens_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
fig = plt.figure(dpi = 200) # I DEFINITELY shouldn't have to run this again
temperature, total_col_density = rotational_diagram(319, 220, plot = True, save = True)
# rotational_diagram(250, 300)
# rotational_diagram(255, 300)
# rotational_diagram(343, 278)
# rotational_diagram(345, 280)
# rotational_diagram(325, 281)
## Step 8: Generate synthetic spectrum for a pixel, and compare to original spectrum

# # Retrieve and plot original spectrum for one of these pixels (takes a while to run)
# test_pixel_x, test_pixel_y = 250, 300
# spectrum = ch3cncube[:,test_pixel_y,test_pixel_x]

# # Plot original spectrum for one of these pixels
# spectrum.quicklook()

# # Print temp/col density inputs from fit for this pixel
# rotational_diagram(test_pixel_x, test_pixel_y)

# %%capture
# # Generate synthetic spectrum for this cube, making sure to use temp/col density inputs from fit
# %run /blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py

# pl.plot(sp.xarr, sp.data)

# pl.plot(spectrum.value-np.array(sp.data))
# pl.title("Residuals (data minus synthetic), 1st attempt")
# pl.xlabel("Spectral units :)")
# pl.ylabel('Data minus model (maybe not same units?)')
from astropy import units as u

fmin = 147.1*u.GHz
fmax = 147.2*u.GHz
import numpy as np

sp_axis = np.linspace(fmin, fmax, 1000)
fillingfactor = 1
offset = 0
species = 'CH3CN'
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule

freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
#temp = 38.32444*u.K
#N_tot = (10**(13.71651))*u.cm**-2
N_tot = total_col_density
total_col_density_map
#[Out]# array([[0., 0., 0., ..., 0., 0., 0.],
#[Out]#        [0., 0., 0., ..., 0., 0., 0.],
#[Out]#        [0., 0., 0., ..., 0., 0., 0.],
#[Out]#        ...,
#[Out]#        [0., 0., 0., ..., 0., 0., 0.],
#[Out]#        [0., 0., 0., ..., 0., 0., 0.],
#[Out]#        [0., 0., 0., ..., 0., 0., 0.]])
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
#temp = 38.32444*u.K
#N_tot = (10**(13.71651))*u.cm**-2
N_tot = total_col_density
temp = temperature
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
pixel_x, pixel_y = 319, 220
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")
plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.spectral_axis, masked_cube[:,pixel_y,pixel_x], 'x')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.spectral_axis, masked_cube.include()[:,pixel_y,pixel_x], 'x')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.spectral_axis, masked_cube.mask.include()[:,pixel_y,pixel_x], 'x')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
masked_cube.spectral_axis
#[Out]# <Quantity [-9343.18661444, -8261.2231362 , -7179.25965795, -6097.29617971,
#[Out]#            -5015.33270146, -3933.36922322, -2851.40574497, -1769.44226673,
#[Out]#             -687.47878848,   394.48468977,  1476.44816801,  2558.41164626,
#[Out]#             3640.3751245 ,  4722.33860275,  5804.30208099,  6886.26555924,
#[Out]#             7968.22903748,  9050.19251573, 10132.15599397, 11214.11947222,
#[Out]#            12296.08295046, 13378.04642871, 14460.00990695, 15541.9733852 ,
#[Out]#            16623.93686344, 17705.90034169, 18787.86381993, 19869.82729818,
#[Out]#            20951.79077642, 22033.75425467, 23115.71773291, 24197.68121116,
#[Out]#            25279.6446894 , 26361.60816765, 27443.5716459 , 28525.53512414,
#[Out]#            29607.49860239] m / s>
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], 'x')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
masked_cube.with_spectral_unit(u.GHz).spectral_axis
#[Out]# <Quantity [135.30202763, 135.30153933, 135.30105104, 135.30056274,
#[Out]#            135.30007444, 135.29958615, 135.29909785, 135.29860956,
#[Out]#            135.29812126, 135.29763297, 135.29714467, 135.29665638,
#[Out]#            135.29616808, 135.29567979, 135.29519149, 135.29470319,
#[Out]#            135.2942149 , 135.2937266 , 135.29323831, 135.29275001,
#[Out]#            135.29226172, 135.29177342, 135.29128513, 135.29079683,
#[Out]#            135.29030854, 135.28982024, 135.28933194, 135.28884365,
#[Out]#            135.28835535, 135.28786706, 135.28737876, 135.28689047,
#[Out]#            135.28640217, 135.28591388, 135.28542558, 135.28493729,
#[Out]#            135.28444899] GHz>
i=0
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')

masked_cube.with_spectral_unit(u.GHz).spectral_axis
#[Out]# <Quantity [135.30202763, 135.30153933, 135.30105104, 135.30056274,
#[Out]#            135.30007444, 135.29958615, 135.29909785, 135.29860956,
#[Out]#            135.29812126, 135.29763297, 135.29714467, 135.29665638,
#[Out]#            135.29616808, 135.29567979, 135.29519149, 135.29470319,
#[Out]#            135.2942149 , 135.2937266 , 135.29323831, 135.29275001,
#[Out]#            135.29226172, 135.29177342, 135.29128513, 135.29079683,
#[Out]#            135.29030854, 135.28982024, 135.28933194, 135.28884365,
#[Out]#            135.28835535, 135.28786706, 135.28737876, 135.28689047,
#[Out]#            135.28640217, 135.28591388, 135.28542558, 135.28493729,
#[Out]#            135.28444899] GHz>
freqs
#[Out]# <Quantity [147103.738 , 147129.2302, 147149.0683, 147163.2441, 147171.7519,
#[Out]#            147174.5883] MHz>
ch3cn_freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
masked_cube
#[Out]# SpectralCube with shape=(37, 512, 512) and unit=K:
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:     37  type_s: VRAD      unit_s: m / s  range:    -9343.187 m / s:   29607.499 m / s
masked_cube.wcs.wcs
#[Out]#        flag: 137
#[Out]#       naxis: 3
#[Out]#       crpix: 0x555a58b2d700
#[Out]#                257.00       257.00      -316.00    
#[Out]#          pc: 0x555a62372f90
#[Out]#     pc[0][]:   1.0000       0.0000       0.0000    
#[Out]#     pc[1][]:   0.0000       1.0000       0.0000    
#[Out]#     pc[2][]:   0.0000       0.0000       1.0000    
#[Out]#       cdelt: 0x555a60d06ba0
#[Out]#               -5.5556e-05   5.5556e-05   1082.0    
#[Out]#       crval: 0x555a5d9a4b20
#[Out]#                266.54      -28.705      -3.5233e+05
#[Out]#       cunit: 0x555a642766b0
#[Out]#              "deg"
#[Out]#              "deg"
#[Out]#              "m/s"
#[Out]#       ctype: 0x555a6243f400
#[Out]#              "RA---SIN"
#[Out]#              "DEC--SIN"
#[Out]#              "VRAD"
#[Out]#     lonpole: 180.000000
#[Out]#     latpole: -28.704931
#[Out]#     restfrq: 135297811000.000000
#[Out]#     restwav: 0.000000
#[Out]#         npv: 0
#[Out]#      npvmax: 64
#[Out]#          pv: 0x555a614964f0
#[Out]#         nps: 0
#[Out]#      npsmax: 8
#[Out]#          ps: 0x555a6148fec0
#[Out]#          cd: 0x555a65237150
#[Out]#     cd[0][]:   0.0000       0.0000       0.0000    
#[Out]#     cd[1][]:   0.0000       0.0000       0.0000    
#[Out]#     cd[2][]:   0.0000       0.0000       0.0000    
#[Out]#       crota: 0x555a5a33bf10
#[Out]#                0.0000       0.0000       0.0000    
#[Out]#      altlin: 0
#[Out]#      velref: 0
#[Out]#         alt: ' '
#[Out]#      colnum: 0
#[Out]#       colax: 0x555a633d37f0
#[Out]#                  0      0      0
#[Out]#       cname: 0x555a64235450
#[Out]#              UNDEFINED
#[Out]#              UNDEFINED
#[Out]#              UNDEFINED
#[Out]#       crder: 0x555a646c1c20
#[Out]#                UNDEFINED    UNDEFINED    UNDEFINED
#[Out]#       csyer: 0x555a6244f080
#[Out]#                UNDEFINED    UNDEFINED    UNDEFINED
#[Out]#       czphs: 0x555a623a7930
#[Out]#                UNDEFINED    UNDEFINED    UNDEFINED
#[Out]#       cperi: 0x555a59d87cc0
#[Out]#                UNDEFINED    UNDEFINED    UNDEFINED
#[Out]#     wcsname: UNDEFINED
#[Out]#     timesys: UNDEFINED
#[Out]#     trefpos: UNDEFINED
#[Out]#     trefdir: UNDEFINED
#[Out]#     plephem: UNDEFINED
#[Out]#    timeunit: UNDEFINED
#[Out]#     dateref: "1858-11-17"
#[Out]#      mjdref:      0.000000000     0.000000000
#[Out]#    timeoffs: UNDEFINED
#[Out]#     dateobs: "2019-11-09T20:04:58.368000"
#[Out]#     datebeg: UNDEFINED
#[Out]#     dateavg: UNDEFINED
#[Out]#     dateend: UNDEFINED
#[Out]#      mjdobs:  58796.836786667
#[Out]#      mjdbeg: UNDEFINED
#[Out]#      mjdavg: UNDEFINED
#[Out]#      mjdend: UNDEFINED
#[Out]#      jepoch: UNDEFINED
#[Out]#      bepoch: UNDEFINED
#[Out]#      tstart: UNDEFINED
#[Out]#       tstop: UNDEFINED
#[Out]#     xposure: UNDEFINED
#[Out]#     telapse: UNDEFINED
#[Out]#     timsyer: UNDEFINED
#[Out]#     timrder: UNDEFINED
#[Out]#     timedel: UNDEFINED
#[Out]#    timepixr: UNDEFINED
#[Out]#      obsgeo:   2225142.180269 -5440307.370349 -2481029.851873
#[Out]#                    -67.754929      -23.022886     5053.796333
#[Out]#    obsorbit: UNDEFINED
#[Out]#     radesys: "ICRS"
#[Out]#     equinox: UNDEFINED
#[Out]#     specsys: "LSRK"
#[Out]#     ssysobs: UNDEFINED
#[Out]#     velosys: UNDEFINED
#[Out]#     zsource: UNDEFINED
#[Out]#     ssyssrc: UNDEFINED
#[Out]#     velangl: UNDEFINED
#[Out]#         aux: 0x0
#[Out]#        ntab: 0
#[Out]#         tab: 0x0
#[Out]#        nwtb: 0
#[Out]#         wtb: 0x0
#[Out]#       types: 0x555a59d8e2f0
#[Out]#             2200 2201 3000
#[Out]#      lngtyp: "RA"
#[Out]#      lattyp: "DEC"
#[Out]#         lng: 0
#[Out]#         lat: 1
#[Out]#        spec: 2
#[Out]#    cubeface: -1
#[Out]#         err: 0x0
#[Out]#         lin: (see below)
#[Out]#         cel: (see below)
#[Out]#         spc: (see below)
#[Out]#      m_flag: 137
#[Out]#     m_naxis: 3
#[Out]#     m_crpix: 0x555a58b2d700  (= crpix)
#[Out]#        m_pc: 0x555a62372f90  (= pc)
#[Out]#     m_cdelt: 0x555a60d06ba0  (= cdelt)
#[Out]#     m_crval: 0x555a5d9a4b20  (= crval)
#[Out]#     m_cunit: 0x555a642766b0  (= cunit)
#[Out]#     m_ctype: 0x555a6243f400  (= ctype)
#[Out]#        m_pv: 0x555a614964f0  (= pv)
#[Out]#        m_ps: 0x555a6148fec0  (= ps)
#[Out]#        m_cd: 0x555a65237150  (= cd)
#[Out]#     m_crota: 0x555a5a33bf10  (= crota)
#[Out]# 
#[Out]#     m_colax: 0x555a633d37f0  (= colax)
#[Out]#     m_cname: 0x555a64235450  (= cname)
#[Out]#     m_crder: 0x555a646c1c20  (= crder)
#[Out]#     m_csyer: 0x555a6244f080  (= csyer)
#[Out]#     m_czphs: 0x555a623a7930  (= czphs)
#[Out]#     m_cperi: 0x555a59d87cc0  (= cperi)
#[Out]#       m_aux: 0x0  (= aux)
#[Out]#       m_tab: 0x0  (= tab)
#[Out]#       m_wtb: 0x0  (= wtb)
#[Out]# 
#[Out]#    lin.*
#[Out]#        flag: 137
#[Out]#       naxis: 3
#[Out]#       crpix: 0x555a58b2d700
#[Out]#                257.00       257.00      -316.00    
#[Out]#          pc: 0x555a62372f90
#[Out]#     pc[0][]:   1.0000       0.0000       0.0000    
#[Out]#     pc[1][]:   0.0000       1.0000       0.0000    
#[Out]#     pc[2][]:   0.0000       0.0000       1.0000    
#[Out]#       cdelt: 0x555a60d06ba0
#[Out]#               -5.5556e-05   5.5556e-05   1082.0    
#[Out]#      dispre: 0x0
#[Out]#      disseq: 0x0
#[Out]#      piximg: (nil)
#[Out]#      imgpix: (nil)
#[Out]#     i_naxis: 0
#[Out]#       unity: 1
#[Out]#      affine: 1
#[Out]#      simple: 1
#[Out]#         err: 0x0
#[Out]#      tmpcrd: 0x555a614fd660
#[Out]#      m_flag: 0
#[Out]#     m_naxis: 0
#[Out]#     m_crpix: 0x0
#[Out]#        m_pc: 0x0
#[Out]#     m_cdelt: 0x0
#[Out]#    m_dispre: 0x0
#[Out]#    m_disseq: 0x0
#[Out]# 
#[Out]#    cel.*
#[Out]#       flag: 137
#[Out]#      offset: 0
#[Out]#        phi0:  0.000000
#[Out]#      theta0: 90.000000
#[Out]#         ref:   266.54      -28.705       180.00      -28.705    
#[Out]#         prj: (see below)
#[Out]#       euler:   266.54       118.70       180.00      -0.48030      0.87710   
#[Out]#     latpreq: 0 (not required)
#[Out]#      isolat: 0
#[Out]#         err: 0x0
#[Out]# 
#[Out]#    prj.*
#[Out]#        flag: 105
#[Out]#        code: "SIN"
#[Out]#          r0: 57.295780
#[Out]#          pv: (0)
#[Out]#               0.0000       0.0000    
#[Out]#        phi0:  0.000000
#[Out]#      theta0: 90.000000
#[Out]#      bounds: 7
#[Out]# 
#[Out]#        name: "orthographic/synthesis"
#[Out]#    category: 1 (zenithal)
#[Out]#     pvrange: 102
#[Out]#   simplezen: 1
#[Out]#   equiareal: 0
#[Out]#   conformal: 0
#[Out]#      global: 0
#[Out]#   divergent: 0
#[Out]#          x0: 0.000000
#[Out]#          y0: 0.000000
#[Out]#         err: 0x0
#[Out]#         w[]:   0.017453     0.0000       1.0000      -1.0000       0.0000    
#[Out]#                0.0000       0.0000       0.0000       0.0000       0.0000    
#[Out]#           m: 0
#[Out]#           n: 0
#[Out]#      prjx2s: 0x2b3cf9d554a0
#[Out]#      prjs2x: 0x2b3cf9d55d50
#[Out]# 
#[Out]#    spc.*
#[Out]#        flag: 0
#[Out]#        type: "    "
#[Out]#        code: "   "
#[Out]#       crval: UNDEFINED
#[Out]#     restfrq: 0.000000
#[Out]#     restwav: 0.000000
#[Out]#          pv: (not used)
#[Out]#           w:   0.0000       0.0000       0.0000      (remainder unused)
#[Out]#     isGrism: 0
#[Out]#         err: 0x0
#[Out]#      spxX2P: 0x0
#[Out]#      spxP2S: 0x0
#[Out]#      spxS2P: 0x0
#[Out]#      spxP2X: 0x0
masked_cube.wcs.wcs.restfrq
#[Out]# 135297811000.0
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], 'x')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-')

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=0.5,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 8):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
ch3cn_freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs7-[i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[7-i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
#temp = 38.32444*u.K
#N_tot = (10**(13.71651))*u.cm**-2
N_tot = total_col_density / 2
temp = temperature
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
pixel_x, pixel_y = 319, 220
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[7-i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
i=0
masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')

masked_cube.with_spectral_unit(u.GHz).spectral_axis
#[Out]# <Quantity [135.30202763, 135.30153933, 135.30105104, 135.30056274,
#[Out]#            135.30007444, 135.29958615, 135.29909785, 135.29860956,
#[Out]#            135.29812126, 135.29763297, 135.29714467, 135.29665638,
#[Out]#            135.29616808, 135.29567979, 135.29519149, 135.29470319,
#[Out]#            135.2942149 , 135.2937266 , 135.29323831, 135.29275001,
#[Out]#            135.29226172, 135.29177342, 135.29128513, 135.29079683,
#[Out]#            135.29030854, 135.28982024, 135.28933194, 135.28884365,
#[Out]#            135.28835535, 135.28786706, 135.28737876, 135.28689047,
#[Out]#            135.28640217, 135.28591388, 135.28542558, 135.28493729,
#[Out]#            135.28444899] GHz>
mod
#[Out]# array([6.25460043e-10, 4.99326042e-10, 3.91325525e-10, 3.01065941e-10,
#[Out]#        2.27377264e-10, 1.68578771e-10, 1.22696283e-10, 8.76641543e-11,
#[Out]#        6.14853377e-11, 4.23342900e-11, 2.86141371e-11, 1.89871636e-11,
#[Out]#        1.23671905e-11, 7.90814150e-12, 4.96212037e-12, 3.06035521e-12,
#[Out]#        1.85221988e-12, 1.09761173e-12, 6.40273480e-13, 3.65870544e-13,
#[Out]#        2.05802172e-13, 1.14334535e-13, 6.09784157e-14, 3.43003573e-14,
#[Out]#        1.52446025e-14, 7.62230093e-15, 3.81115029e-15, 3.81115012e-15,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 3.81111997e-15, 3.81111980e-15,
#[Out]#        1.14333589e-14, 2.28667167e-14, 4.57334314e-14, 9.14668586e-14,
#[Out]#        1.79122590e-13, 3.35378451e-13, 6.25023449e-13, 1.13571329e-12,
#[Out]#        2.03513714e-12, 3.57482875e-12, 6.16638877e-12, 1.04386515e-11,
#[Out]#        1.73520182e-11, 2.83127911e-11, 4.53522948e-11, 7.13174332e-11,
#[Out]#        1.10091733e-10, 1.66831633e-10, 2.48187530e-10, 3.62448599e-10,
#[Out]#        5.19619013e-10, 7.31292184e-10, 1.01034586e-09, 1.37030569e-09,
#[Out]#        1.82446864e-09, 2.38464534e-09, 3.05972331e-09, 3.85398621e-09,
#[Out]#        4.76549798e-09, 5.78464317e-09, 6.89310174e-09, 8.06349082e-09,
#[Out]#        9.25979155e-09, 1.04387745e-08, 1.15522824e-08, 1.25503553e-08,
#[Out]#        1.33848781e-08, 1.40133758e-08, 1.44026079e-08, 1.45314686e-08,
#[Out]#        1.43928883e-08, 1.39944663e-08, 1.33577965e-08, 1.25165119e-08,
#[Out]#        1.15133542e-08, 1.03965760e-08, 9.21613703e-09, 8.02005737e-09,
#[Out]#        6.85134922e-09, 5.74572269e-09, 4.73024047e-09, 3.82289098e-09,
#[Out]#        3.03298858e-09, 2.36221394e-09, 1.80608518e-09, 1.35558520e-09,
#[Out]#        9.98815703e-10, 7.22460653e-10, 5.12994412e-10, 3.57588779e-10,
#[Out]#        2.44692272e-10, 1.64373136e-10, 1.08395561e-10, 7.01701347e-11,
#[Out]#        4.45937817e-11, 2.78210915e-11, 1.70394650e-11, 1.02442585e-11,
#[Out]#        6.04441715e-12, 3.50240801e-12, 1.99320926e-12, 1.11284336e-12,
#[Out]#        6.09777155e-13, 3.27755206e-13, 1.71499809e-13, 9.14665609e-14,
#[Out]#        4.57332784e-14, 2.28666382e-14, 1.14333186e-14, 3.81110602e-15,
#[Out]#        3.81110584e-15, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 3.81108620e-15, 7.62217206e-15,
#[Out]#        1.90554293e-14, 4.19219425e-14, 8.76549668e-14, 1.75309926e-13,
#[Out]#        3.50619836e-13, 6.93617469e-13, 1.34150186e-12, 2.54580455e-12,
#[Out]#        4.74098909e-12, 8.67402788e-12, 1.55682787e-11, 2.74398045e-11,
#[Out]#        4.74784819e-11, 8.06425294e-11, 1.34458835e-10, 2.20082437e-10,
#[Out]#        3.53641834e-10, 5.57843453e-10, 8.63839073e-10, 1.31318092e-09,
#[Out]#        1.95969669e-09, 2.87094923e-09, 4.12889582e-09, 5.82928223e-09,
#[Out]#        8.07919599e-09, 1.09924207e-08, 1.46821910e-08, 1.92513223e-08,
#[Out]#        2.47800643e-08, 3.13124082e-08, 3.88420999e-08, 4.73000025e-08,
#[Out]#        5.65446828e-08, 6.63581964e-08, 7.64486070e-08, 8.64603026e-08,
#[Out]#        9.59922396e-08, 1.04623126e-07, 1.11941583e-07, 1.17578354e-07,
#[Out]#        1.21237087e-07, 1.22720133e-07, 1.21946212e-07, 1.18957821e-07,
#[Out]#        1.13917352e-07, 1.07092472e-07, 9.88325926e-08, 8.95392864e-08,
#[Out]#        7.96341354e-08, 6.95275738e-08, 5.95918788e-08, 5.01405727e-08,
#[Out]#        4.14155710e-08, 3.35822858e-08, 2.67318475e-08, 2.08891080e-08,
#[Out]#        1.60244474e-08, 1.20675284e-08, 8.92125222e-09, 6.47449167e-09,
#[Out]#        4.61272382e-09, 3.22612864e-09, 2.21501976e-09, 1.49295433e-09,
#[Out]#        9.87842029e-10, 6.41651585e-10, 4.09153146e-10, 2.56119438e-10,
#[Out]#        1.57389740e-10, 9.49452858e-11, 5.62285828e-11, 3.26875785e-11,
#[Out]#        1.86552046e-11, 1.04537740e-11, 5.74709825e-12, 3.10221337e-12,
#[Out]#        1.64638343e-12, 8.57491333e-13, 4.38273328e-13, 2.21042190e-13,
#[Out]#        1.06710018e-13, 5.33550066e-14, 2.28664304e-14, 1.14332147e-14,
#[Out]#        3.81107139e-15, 3.81107121e-15, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 3.81106191e-15,
#[Out]#        7.62212347e-15, 1.52442463e-14, 3.81106139e-14, 7.62212244e-14,
#[Out]#        1.56253503e-13, 3.20129113e-13, 6.36447137e-13, 1.24621679e-12,
#[Out]#        2.39334590e-12, 4.50848420e-12, 8.34241037e-12, 1.51451518e-11,
#[Out]#        2.69899246e-11, 4.72228382e-11, 8.11107756e-11, 1.36763669e-10,
#[Out]#        2.26380714e-10, 3.67862451e-10, 5.86819186e-10, 9.18956716e-10,
#[Out]#        1.41273644e-09, 2.13206216e-09, 3.15873440e-09, 4.59408914e-09,
#[Out]#        6.55933700e-09, 9.19377589e-09, 1.26503244e-08, 1.70876786e-08,
#[Out]#        2.26588713e-08, 2.94962796e-08, 3.76938038e-08, 4.72875001e-08,
#[Out]#        5.82366701e-08, 7.04077437e-08, 8.35637722e-08, 9.73619814e-08,
#[Out]#        1.11361352e-07, 1.25041256e-07, 1.37830679e-07, 1.49146198e-07,
#[Out]#        1.58435398e-07, 1.65221282e-07, 1.69142792e-07, 1.69986617e-07,
#[Out]#        1.67706425e-07, 1.62427067e-07, 1.54433264e-07, 1.44144152e-07,
#[Out]#        1.32076921e-07, 1.18803866e-07, 1.04907843e-07, 9.09408623e-08,
#[Out]#        7.73898330e-08, 6.46520798e-08, 5.30218424e-08, 4.26875163e-08,
#[Out]#        3.37381046e-08, 2.61766560e-08, 1.99379936e-08, 1.49081058e-08,
#[Out]#        1.09430187e-08, 7.88543466e-09, 5.57812110e-09, 3.87368169e-09,
#[Out]#        2.64078756e-09, 1.76732514e-09, 1.16111271e-09, 7.48867552e-10,
#[Out]#        4.74144166e-10, 2.94704674e-10, 1.79820568e-10, 1.07711692e-10,
#[Out]#        6.33396399e-11, 3.65632056e-11, 2.07206733e-11, 1.15284227e-11,
#[Out]#        6.29585235e-12, 3.37658894e-12, 1.77594851e-12, 9.18462598e-13,
#[Out]#        4.64947850e-13, 2.32473915e-13, 1.14331428e-13, 5.71657115e-14,
#[Out]#        3.04883781e-14, 2.66773296e-14, 3.81104692e-14, 6.85988414e-14,
#[Out]#        1.41008723e-13, 2.93450573e-13, 5.94523212e-13, 1.18142428e-12,
#[Out]#        2.30949381e-12, 4.43224616e-12, 8.34618973e-12, 1.54309227e-11,
#[Out]#        2.80073711e-11, 4.99018235e-11, 8.72919823e-11, 1.49896009e-10,
#[Out]#        2.52683684e-10, 4.18159217e-10, 6.79333866e-10, 1.08342264e-09,
#[Out]#        1.69624225e-09, 2.60706636e-09, 3.93360663e-09, 5.82645638e-09,
#[Out]#        8.47213180e-09, 1.20935979e-08, 1.69470098e-08, 2.33133824e-08,
#[Out]#        3.14841611e-08, 4.17401203e-08, 5.43237754e-08, 6.94066269e-08,
#[Out]#        8.70535916e-08, 1.07188277e-07, 1.29563563e-07, 1.53742407e-07,
#[Out]#        1.79093594e-07, 2.04805978e-07, 2.29923025e-07, 2.53396932e-07,
#[Out]#        2.74158877e-07, 2.91199379e-07, 3.03650742e-07, 3.10862649e-07,
#[Out]#        3.12462220e-07, 3.08391604e-07, 2.98918779e-07, 2.84620722e-07,
#[Out]#        2.66341575e-07, 2.45131482e-07, 2.22173801e-07, 1.98709067e-07,
#[Out]#        1.75963523e-07, 1.55088192e-07, 1.37111832e-07, 1.22908419e-07,
#[Out]#        1.13177154e-07, 1.08431405e-07, 1.08992240e-07, 1.14982762e-07,
#[Out]#        1.26320958e-07, 1.42710830e-07, 1.63634062e-07, 1.88346168e-07,
#[Out]#        2.15882405e-07, 2.45078564e-07, 2.74610472e-07, 3.03053530e-07,
#[Out]#        3.28960516e-07, 3.50952333e-07, 3.67813514e-07, 3.78582408e-07,
#[Out]#        3.82625460e-07, 3.79686316e-07, 3.69903248e-07, 3.53792255e-07,
#[Out]#        3.32197616e-07, 3.06215659e-07, 2.77100833e-07, 2.46164797e-07,
#[Out]#        2.14679533e-07, 1.83793971e-07, 1.54471162e-07, 1.27449715e-07,
#[Out]#        1.03229899e-07, 8.20818680e-08, 6.40713671e-08, 4.90970905e-08,
#[Out]#        3.69336763e-08, 2.72749686e-08, 1.97733858e-08, 1.40725466e-08,
#[Out]#        9.83194274e-09, 6.74343303e-09, 4.54043573e-09, 3.00116015e-09,
#[Out]#        1.94740240e-09, 1.24050188e-09, 7.75735233e-10, 4.76218712e-10,
#[Out]#        2.86993428e-10, 1.69789002e-10, 9.86103931e-11, 5.62241184e-11,
#[Out]#        3.14676693e-11, 1.72906393e-11, 9.32558796e-12, 4.93909336e-12,
#[Out]#        2.56863331e-12, 1.31099380e-12, 6.55496869e-13, 3.23937392e-13,
#[Out]#        1.56252146e-13, 7.24095280e-14, 3.42992486e-14, 1.52441098e-14,
#[Out]#        7.62205455e-15, 3.81102710e-15, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
#temp = 38.32444*u.K
#N_tot = (10**(13.71651))*u.cm**-2
N_tot = total_col_density - np.log10(0.5)
temp = temperature
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
pixel_x, pixel_y = 319, 220
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[7-i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
v_cen = 10*u.km/u.s
v_disp = 1.5*u.km/u.s
#temp = 38.32444*u.K
#N_tot = (10**(13.71651))*u.cm**-2
N_tot = total_col_density - np.log10(2)
temp = temperature
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
pixel_x, pixel_y = 319, 220
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150
# fig = plt.figure(dpi = display_dpi)
mod_sp.plotter(label = 'Synthetic spectrum', color = 'tab:red', linewidth = 1)
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, ymin = min(mod) - 0.3, ymax = max(mod) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

for i in range(0, 5):

    # Import masked cube and get channel width
    masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{i}_masked.fits', format='fits')
    plt.plot(masked_cube.with_spectral_unit(u.GHz, rest_value=ch3cn_freqs[7-i]).spectral_axis,
             masked_cube.mask.include()[:,pixel_y,pixel_x], '-', alpha=0.5, linewidth=1,
            label=str(i))

plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend() # Is this hacking???
data_sp_K_pyspeckit.plotter.savefig('../figures/ch3cn_spec_datasynth.pdf')
mod
#[Out]# array([2.25635208e-03, 1.80134176e-03, 1.41173178e-03, 1.08611322e-03,
#[Out]#        8.20285192e-04, 6.08165168e-04, 4.42634396e-04, 3.16253810e-04,
#[Out]#        2.21816277e-04, 1.52727815e-04, 1.03230963e-04, 6.84965922e-05,
#[Out]#        4.46164615e-05, 2.85291264e-05, 1.79080800e-05, 1.10351133e-05,
#[Out]#        6.67531540e-06, 3.96400336e-06, 2.31080591e-06, 1.32239184e-06,
#[Out]#        7.42889278e-07, 4.09689899e-07, 2.21795976e-07, 1.17874348e-07,
#[Out]#        6.14967533e-08, 3.14957781e-08, 1.58350589e-08, 7.81545695e-09,
#[Out]#        3.78665950e-09, 1.80104648e-09, 8.40933971e-10, 3.85448220e-10,
#[Out]#        1.73433969e-10, 7.66079078e-11, 3.32179739e-11, 1.41393618e-11,
#[Out]#        5.90728028e-12, 2.42389038e-12, 9.75653945e-13, 3.84925953e-13,
#[Out]#        1.48634767e-13, 5.71672156e-14, 2.28668852e-14, 7.62229472e-15,
#[Out]#        3.81114719e-15, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 3.81112273e-15, 1.52444902e-14,
#[Out]#        3.81112238e-14, 1.02900300e-13, 2.82023031e-13, 7.43168764e-13,
#[Out]#        1.92461646e-12, 4.88966891e-12, 1.22108328e-11, 2.99249235e-11,
#[Out]#        7.19958869e-11, 1.70044589e-10, 3.94256621e-10, 8.97354997e-10,
#[Out]#        2.00503421e-09, 4.39794118e-09, 9.46993569e-09, 2.00177124e-08,
#[Out]#        4.15386086e-08, 8.46173591e-08, 1.69214279e-07, 3.32188480e-07,
#[Out]#        6.40180070e-07, 1.21112719e-06, 2.24929999e-06, 4.10086077e-06,
#[Out]#        7.33960490e-06, 1.28955642e-05, 2.22422139e-05, 3.76604637e-05,
#[Out]#        6.25984061e-05, 1.02143520e-04, 1.63616933e-04, 2.57285608e-04,
#[Out]#        3.97166374e-04, 6.01864964e-04, 8.95354694e-04, 1.30755643e-03,
#[Out]#        1.87454041e-03, 2.63814070e-03, 3.64476547e-03, 4.94321383e-03,
#[Out]#        6.58138172e-03, 8.60185996e-03, 1.10365937e-02, 1.39009681e-02,
#[Out]#        1.71878830e-02, 2.08625454e-02, 2.48588015e-02, 2.90778097e-02,
#[Out]#        3.33897083e-02, 3.76386449e-02, 4.16511360e-02, 4.52472698e-02,
#[Out]#        4.82538160e-02, 5.05179445e-02, 5.19200571e-02, 5.23842360e-02,
#[Out]#        5.18850444e-02, 5.04498318e-02, 4.81562557e-02, 4.51253301e-02,
#[Out]#        4.15108655e-02, 3.74865784e-02, 3.32323682e-02, 2.89212568e-02,
#[Out]#        2.47082752e-02, 2.07222201e-02, 1.70607493e-02, 1.37888319e-02,
#[Out]#        1.09401732e-02, 8.52095399e-03, 6.51507847e-03, 4.89011074e-03,
#[Out]#        3.60317721e-03, 2.60627799e-03, 1.85065025e-03, 1.29002091e-03,
#[Out]#        8.82750974e-04, 5.92992143e-04, 3.91047165e-04, 2.53150600e-04,
#[Out]#        1.60878688e-04, 1.00366295e-04, 6.14677270e-05, 3.69552673e-05,
#[Out]#        2.18109964e-05, 1.26370194e-05, 7.18759825e-06, 4.01321980e-06,
#[Out]#        2.19974385e-06, 1.18364453e-06, 6.25230998e-07, 3.24212479e-07,
#[Out]#        1.65039907e-07, 8.24742243e-08, 4.04592216e-08, 1.94843938e-08,
#[Out]#        9.21141234e-09, 4.27499345e-09, 1.94766928e-09, 8.71092911e-10,
#[Out]#        3.82459647e-10, 1.64845535e-10, 6.97470291e-11, 2.89720175e-11,
#[Out]#        1.18144238e-11, 4.72958043e-12, 1.85981881e-12, 7.16487542e-13,
#[Out]#        2.70588368e-13, 9.90886937e-14, 3.81110343e-14, 1.14333098e-14,
#[Out]#        3.81110309e-15, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 3.81108896e-15, 1.14332664e-14,
#[Out]#        3.42997975e-14, 9.90882995e-14, 2.70587267e-13, 7.46973267e-13,
#[Out]#        2.01225442e-12, 5.32790068e-12, 1.38494923e-11, 3.53211581e-11,
#[Out]#        8.84515237e-11, 2.17426328e-10, 5.24683766e-10, 1.24296307e-09,
#[Out]#        2.89061768e-09, 6.59925049e-09, 1.47900938e-08, 3.25401311e-08,
#[Out]#        7.02813205e-08, 1.49015928e-07, 3.10168487e-07, 6.33774605e-07,
#[Out]#        1.27128869e-06, 2.50337388e-06, 4.83926531e-06, 9.18343852e-06,
#[Out]#        1.71081630e-05, 3.12876960e-05, 5.61714852e-05, 9.89988857e-05,
#[Out]#        1.71283965e-04, 2.90920950e-04, 4.85070656e-04, 7.93974348e-04,
#[Out]#        1.27578953e-03, 2.01243644e-03, 3.11627070e-03, 4.73714956e-03,
#[Out]#        7.06914199e-03, 1.03557715e-02, 1.48923246e-02, 2.10234806e-02,
#[Out]#        2.91344212e-02, 3.96337555e-02, 5.29271527e-02, 6.93815452e-02,
#[Out]#        8.92811367e-02, 1.12778097e-01, 1.39842531e-01, 1.70217745e-01,
#[Out]#        2.03387693e-01, 2.38563371e-01, 2.74693745e-01, 3.10504430e-01,
#[Out]#        3.44564071e-01, 3.75374639e-01, 4.01478208e-01, 4.21569840e-01,
#[Out]#        4.34604616e-01, 4.39886744e-01, 4.37130389e-01, 4.26484982e-01,
#[Out]#        4.08521967e-01, 3.84184580e-01, 3.54706660e-01, 3.21510024e-01,
#[Out]#        2.86092096e-01, 2.49915976e-01, 2.14313986e-01, 1.80413150e-01,
#[Out]#        1.49087601e-01, 1.20939103e-01, 9.63034298e-02, 7.52776378e-02,
#[Out]#        5.77617175e-02, 4.35076573e-02, 3.21695509e-02, 2.33496581e-02,
#[Out]#        1.66369954e-02, 1.16367207e-02, 7.99005827e-03, 5.38560968e-03,
#[Out]#        3.56358432e-03, 2.31476958e-03, 1.47603986e-03, 9.23971937e-04,
#[Out]#        5.67793510e-04, 3.42525909e-04, 2.02846761e-04, 1.17927358e-04,
#[Out]#        6.73027845e-05, 3.77071363e-05, 2.07389219e-05, 1.11974967e-05,
#[Out]#        5.93509770e-06, 3.08821106e-06, 1.57745963e-06, 7.91009535e-07,
#[Out]#        3.89383343e-07, 1.88167757e-07, 8.92658261e-08, 4.15716678e-08,
#[Out]#        1.90056110e-08, 8.52979096e-09, 3.75808953e-09, 1.62542935e-09,
#[Out]#        6.90146792e-10, 2.87663414e-10, 1.17708719e-10, 4.72801366e-11,
#[Out]#        1.86437545e-11, 7.21816627e-12, 2.74397016e-12, 1.02136662e-12,
#[Out]#        3.73484793e-13, 1.33387420e-13, 4.57328277e-14, 1.52442752e-14,
#[Out]#        3.81106863e-15, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 3.81106467e-15,
#[Out]#        7.62212899e-15, 2.66774502e-14, 7.24102188e-14, 2.01986391e-13,
#[Out]#        5.60226379e-13, 1.53204758e-12, 4.10070428e-12, 1.07738759e-11,
#[Out]#        2.77940833e-11, 7.03712772e-11, 1.74923970e-10, 4.26846633e-10,
#[Out]#        1.02250805e-09, 2.40455924e-09, 5.55109775e-09, 1.25804411e-08,
#[Out]#        2.79889481e-08, 6.11295151e-08, 1.31065702e-07, 2.75867737e-07,
#[Out]#        5.70015275e-07, 1.15623421e-06, 2.30238986e-06, 4.50075779e-06,
#[Out]#        8.63706408e-06, 1.62712289e-05, 3.00918066e-05, 5.46323444e-05,
#[Out]#        9.73699690e-05, 1.70362408e-04, 2.92614558e-04, 4.93391119e-04,
#[Out]#        8.16694604e-04, 1.32709067e-03, 2.11696399e-03, 3.31510659e-03,
#[Out]#        5.09625993e-03, 7.69084805e-03, 1.13936694e-02, 1.65698052e-02,
#[Out]#        2.36555384e-02, 3.31517714e-02, 4.56074320e-02, 6.15908025e-02,
#[Out]#        8.16477052e-02, 1.06247071e-01, 1.35716511e-01, 1.70172873e-01,
#[Out]#        2.09455053e-01, 2.53068031e-01, 3.00147740e-01, 3.49455631e-01,
#[Out]#        3.99409358e-01, 4.48152184e-01, 4.93658798e-01, 5.33870093e-01,
#[Out]#        5.66844849e-01, 5.90913005e-01, 6.04813975e-01, 6.07804421e-01,
#[Out]#        5.99723006e-01, 5.81004571e-01, 5.52642033e-01, 5.16100508e-01,
#[Out]#        4.73193664e-01, 4.25936453e-01, 3.76390535e-01, 3.26518661e-01,
#[Out]#        2.78062007e-01, 2.32450544e-01, 1.90751530e-01, 1.53656072e-01,
#[Out]#        1.21499255e-01, 9.43060618e-02, 7.18537076e-02, 5.37408990e-02,
#[Out]#        3.94557334e-02, 2.84359896e-02, 2.01179253e-02, 1.39719819e-02,
#[Out]#        9.52567616e-03, 6.37527824e-03, 4.18861179e-03, 2.70153368e-03,
#[Out]#        1.71049283e-03, 1.06317239e-03, 6.48721634e-04, 3.88584821e-04,
#[Out]#        2.28500095e-04, 1.31904701e-04, 7.47493909e-05, 4.15842310e-05,
#[Out]#        2.27103369e-05, 1.21756825e-05, 6.40829241e-06, 3.31124859e-06,
#[Out]#        1.68013263e-06, 8.38082422e-07, 4.13195618e-07, 2.06434077e-07,
#[Out]#        1.15732111e-07, 9.46390350e-08, 1.32746752e-07, 2.51608732e-07,
#[Out]#        5.15036987e-07, 1.05772964e-06, 2.14234961e-06, 4.26382401e-06,
#[Out]#        8.33239702e-06, 1.59857639e-05, 3.01075365e-05, 5.56663085e-05,
#[Out]#        1.01037903e-04, 1.80032415e-04, 3.14913593e-04, 5.40762037e-04,
#[Out]#        9.11579878e-04, 1.50853884e-03, 2.45070334e-03, 3.90837646e-03,
#[Out]#        6.11887998e-03, 9.40405896e-03, 1.41880898e-02, 2.10132997e-02,
#[Out]#        3.05507636e-02, 4.36015861e-02, 6.10842255e-02, 8.40032369e-02,
#[Out]#        1.13395679e-01, 1.50253346e-01, 1.95422006e-01, 2.49482800e-01,
#[Out]#        3.12625349e-01, 3.84526325e-01, 4.64250190e-01, 5.50189703e-01,
#[Out]#        6.40061881e-01, 7.30970273e-01, 8.19536926e-01, 9.02098528e-01,
#[Out]#        9.74952162e-01, 1.03462859e+00, 1.07816617e+00, 1.10335735e+00,
#[Out]#        1.10894207e+00, 1.09472814e+00, 1.06162704e+00, 1.01160257e+00,
#[Out]#        9.47539911e-01, 8.73050806e-01, 7.92236882e-01, 7.09436323e-01,
#[Out]#        6.28978458e-01, 5.54966639e-01, 4.91102648e-01, 4.40557273e-01,
#[Out]#        4.05883255e-01, 3.88960544e-01, 3.90960847e-01, 4.12319608e-01,
#[Out]#        4.52708282e-01, 5.11006952e-01, 5.85284973e-01, 6.72803467e-01,
#[Out]#        7.70056308e-01, 8.72864746e-01, 9.76535020e-01, 1.07607923e+00,
#[Out]#        1.16648918e+00, 1.24304306e+00, 1.30161740e+00, 1.33897336e+00,
#[Out]#        1.35298730e+00, 1.34280030e+00, 1.30886973e+00, 1.25291597e+00,
#[Out]#        1.17776867e+00, 1.08712758e+00, 9.85262599e-01, 8.76683607e-01,
#[Out]#        7.65813861e-01, 6.56698821e-01, 5.52776534e-01, 4.56726495e-01,
#[Out]#        3.70402963e-01, 2.94847885e-01, 2.30369658e-01, 1.76668226e-01,
#[Out]#        1.32984981e-01, 9.82572609e-02, 7.12610991e-02, 5.07311015e-02,
#[Out]#        3.54517548e-02, 2.43192329e-02, 1.63763374e-02, 1.08253890e-02,
#[Out]#        7.02480278e-03, 4.47498743e-03, 2.79845900e-03, 1.71797674e-03,
#[Out]#        1.03535209e-03, 6.12536789e-04, 3.55754345e-04, 2.02834844e-04,
#[Out]#        1.13529698e-04, 6.23808276e-05, 3.36486474e-05, 1.78180011e-05,
#[Out]#        9.26243723e-06, 4.72679231e-06, 2.36800558e-06, 1.16459222e-06,
#[Out]#        5.62263665e-07, 2.66490127e-07, 1.23993007e-07, 5.66354350e-08,
#[Out]#        2.53953519e-08, 1.11788020e-08, 4.83071292e-09, 2.04928436e-09,
#[Out]#        8.53429860e-10, 3.48907090e-10, 1.40032348e-10, 5.51722244e-11,
#[Out]#        2.13379340e-11, 8.10224069e-12, 3.01833224e-12, 1.10519736e-12,
#[Out]#        3.96346621e-13, 1.41007926e-13, 4.95433232e-14, 1.52440988e-14,
#[Out]#        3.81102452e-15, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
