########################################################
# Started Logging At: 2022-03-21 17:13:11
########################################################
########################################################
# Started Logging At: 2022-03-21 17:13:11
########################################################
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:13:13
# # Started Logging At: 2022-03-21 17:13:13
########################################################
########################################################
########################################################
# Started Logging At: 2022-03-21 17:13:18
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:13:19
########################################################
########################################################
# Started Logging At: 2022-03-21 17:13:23
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:13:23
########################################################
# Front matter
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150 # For making plots easier to read
plt.rcParams['figure.facecolor'] = 'w' # For making plot axis labels visible even on dark backgrounds
import warnings
warnings.filterwarnings('ignore')
from dask.diagnostics import ProgressBar
ProgressBar().register() # For showing time elapsed on cells
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these... but I don't want other lines creeping in
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
########################################################
# Started Logging At: 2022-03-21 17:13:30
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:13:30
########################################################
get_ipython().run_line_magic('matplotlib', 'inline')
def retrieve_cube(results, freq_spw):
    '''Get methyl cyanide (target molecule) cube'''
    fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
    ch3cncube = SpectralCube.read(fn, format='casa_image')
    ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN',
                                                                                        fmin=fmin, 
                                                                                        fmax=fmax, 
                                                                                        catalog='JPL')
    # We're readying the partition function for use with temperature map later!
    ch3cn_A_new = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij), this line gets us just A_ij and gives it units
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    # ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A_new, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A_orig = retrieve_cube(results, freq_spw)
ch3cncube
from astropy import units as u

fmin = 147.0*u.GHz
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
########################################################
# Started Logging At: 2022-03-21 17:13:37
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:13:39
########################################################
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
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
from spectral_cube import OneDSpectrum
from astropy.io import fits
import pyspeckit

spectra = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/spectra/'
max_fn = spectra+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
mean_fn = spectra+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
# Plot methyl cyanide max spectrum
fig1 = plt.figure(1, figsize = (9, 5))
max_kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn)).to(u.K)
max_kspectrum_ps = pyspeckit.Spectrum.from_hdu(max_kspectrum.hdu)
max_kspectrum_ps.xarr.convert_to_unit('GHz')
max_kspectrum_ps.plotter(figure = fig1, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                         xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                         ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                         title = "CH$_3$CN Max Spectrum")
max_kspectrum_ps.plotter.savefig(f'../figures/ch3cn_max_spec.pdf', facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
v_cen = 50*u.km/u.s
v_disp = 1.5*u.km/u.s
temp = 100*u.K
N_tot = (10**(13.69))*u.cm**-2
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp.unit
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
# Plot methyl cyanide mean spectrum
fig2 = plt.figure(2, figsize = (9, 5))
mean_kspectrum = OneDSpectrum.from_hdu(fits.open(mean_fn)).to(u.K)
mean_kspectrum_ps = pyspeckit.Spectrum.from_hdu(mean_kspectrum.hdu)
mean_kspectrum_ps.xarr.convert_to_unit('GHz')
mean_kspectrum_ps.plotter(figure = fig2, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                          xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                          ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                          title = "CH$_3$CN Mean Spectrum")
plt.show()
freq_spw_2 = '135_spw47'
fmin = 135.25*u.GHz
fmax = 135.32*u.GHz
max_fn_2 = spectra+'source_ab_'+freq_spw_2+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
mean_fn_2 = spectra+'source_ab_'+freq_spw_2+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
# Plot formaldehyde max spectrum
fig3 = plt.figure(3, figsize = (9, 5))
max_kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn_2)).to(u.K)
max_kspectrum_ps = pyspeckit.Spectrum.from_hdu(max_kspectrum.hdu)
max_kspectrum_ps.xarr.convert_to_unit('GHz')
max_kspectrum_ps.plotter(figure = fig1, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                         xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                         ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                         title = "H$_2$CS Max Spectrum")
max_kspectrum_ps.plotter.savefig(f'../figures/h2cs_max_spec.pdf', facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
fn = results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
ch3cncube = SpectralCube.read(fn, format='casa_image')
from pylab import imshow
from astropy.io import fits
from astropy.visualization import quantity_support
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
noise_map = generate_noise_map()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []
    
    kmax = len(ch3cn_freqs) - 1 # Maximum rotational quantum number for the ladder

    for i in range(len(ch3cn_freqs)): # i is the index of the loop
        
        kk = kmax - i # Get rotational quantum number of this iteration

        if kk == 0 or kk == 1:
            # Import k = 0 and k = 1 cubes
            if kk == 0:
                # We need something to be called masked_cube to get channel_width and number_of_pixels
                masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
            elif kk == 1:
                # We need something to be called masked_cube to get channel_width and number_of_pixels
                masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_1_masked.fits', format='fits')
            masked_cube_0 = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
            print("k = 0:", masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]))
            masked_cube_1 = SpectralCube.read(f'methyl_cyanide/ch3cn_1_masked.fits', format='fits')
            print("k = 1:", masked_cube_1.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 1]))
            # Get frequency range we want: highest freq in k = 1, highest freq in k = 0
            highest_freq_1 = np.max(masked_cube_1.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 1]).spectral_axis)
            highest_freq_0 = np.max(masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]).spectral_axis)
            # Get extra part of k = 0 that's not in k = 1 cube
            extra_part_of_0 = masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]).spectral_slab(highest_freq_1, 
                                                                                                                        highest_freq_0)
            print("section of k = 0 that's not in k = 1:", extra_part_of_0)
            extra_part_of_0 = extra_part_of_0.with_spectral_unit(u.km/u.s, velocity_convention = 'radio')
            mom0_sum = extra_part_of_0.moment0() + masked_cube_1.moment0()
            # Split equally between k = 0 and k = 1 component
            mom0 = mom0_sum/2
        else: 
            # Import masked cube and get channel width
            masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{kk}_masked.fits', format='fits')
            # primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb',
            #                                  format='casa_image')
            # masked_cube = masked_cube/primary_beam # Correct for effect of primary beam

            # Calculate moment 0 and moment 1 maps of cube
            mom0 = masked_cube.moment0()
            # mom1 = masked_cube.moment1()
    
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[kmax - kk], ch3cn_A[kmax - kk])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[kmax - kk]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[kmax - kk], ch3cn_A[kmax - kk])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)
        # As of Feb. 23, 2022, I am once again convinced this is correct.

        # Append upper state column density maps and error maps into lists
        log_N_upper_gs.append(log_N_upper_g)
        log_N_upper_g_errs.append(log_N_upper_g_err)

    log_N_upper_gs = np.array(log_N_upper_gs)
    log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    return log_N_upper_gs, log_N_upper_g_errs
log_N_upper_gs, log_N_upper_g_errs = generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g)
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
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s, for treatment as upper limits
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
# # Set the errors for NaN values to be huge so these are ignored in the fit
# ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # Probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # Use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0
# ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # Set errors to be big so that they get ignored in the fit
# # This is my attempt at a brute-force solution to add upper limits. It's way too slow.
# flag = False
# for i in range(len(ch3cn_freqs)):
#     for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
#         while flag == False:
#             if ln_N_upper_gs[i,index[0],index[1]] == 0. and not np.isnan(noise_map[index[0],index[1]].value):
#                 ln_N_upper_g_errs[i,index[0],index[1]] = noise_map[index[0],index[1]].value
#                 flag = True
from scipy.optimize import curve_fit
import tqdm
print(tqdm.__version__) # As of Feb. 28, 2022: 4.62.3
from tqdm.notebook import tqdm
########################################################
# Started Logging At: 2022-03-21 17:14:58
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:14:58
########################################################
# Front matter
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150 # For making plots easier to read
plt.rcParams['figure.facecolor'] = 'w' # For making plot axis labels visible even on dark backgrounds
import warnings
warnings.filterwarnings('ignore')
from dask.diagnostics import ProgressBar
ProgressBar().register() # For showing time elapsed on cells
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these... but I don't want other lines creeping in
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
    ch3cn_A_new = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij), this line gets us just A_ij and gives it units
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    # ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A_new, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A_orig = retrieve_cube(results, freq_spw)
ch3cncube
########################################################
# Started Logging At: 2022-03-21 17:15:06
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:15:07
########################################################
########################################################
# Started Logging At: 2022-03-21 17:15:09
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:15:10
########################################################
# Front matter
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 150 # For making plots easier to read
plt.rcParams['figure.facecolor'] = 'w' # For making plot axis labels visible even on dark backgrounds
import warnings
warnings.filterwarnings('ignore')
from dask.diagnostics import ProgressBar
ProgressBar().register() # For showing time elapsed on cells
from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants

# User inputs
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these... but I don't want other lines creeping in
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
    ch3cn_A_new = 10**ch3cn_A*u.s**-1 # Original is log_10(A_ij), this line gets us just A_ij and gives it units
    # DO NOT RUN THE FOLLOWING LINE, WHATEVER YOU DO
    # ch3cn_E_U = ch3cn_E_U/constants.k_B # Original is in erg
    return ch3cncube, ch3cn_freqs, ch3cn_A_new, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A
ch3cncube, ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc, ch3cn_A_orig = retrieve_cube(results, freq_spw)
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
from spectral_cube import OneDSpectrum
from astropy.io import fits
import pyspeckit

spectra = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/spectra/'
max_fn = spectra+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
mean_fn = spectra+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
# Plot methyl cyanide max spectrum
fig1 = plt.figure(1, figsize = (9, 5))
max_kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn)).to(u.K)
max_kspectrum_ps = pyspeckit.Spectrum.from_hdu(max_kspectrum.hdu)
max_kspectrum_ps.xarr.convert_to_unit('GHz')
max_kspectrum_ps.plotter(figure = fig1, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                         xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                         ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                         title = "CH$_3$CN Max Spectrum")
max_kspectrum_ps.plotter.savefig(f'../figures/ch3cn_max_spec.pdf', facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
# Plot methyl cyanide mean spectrum
fig2 = plt.figure(2, figsize = (9, 5))
mean_kspectrum = OneDSpectrum.from_hdu(fits.open(mean_fn)).to(u.K)
mean_kspectrum_ps = pyspeckit.Spectrum.from_hdu(mean_kspectrum.hdu)
mean_kspectrum_ps.xarr.convert_to_unit('GHz')
mean_kspectrum_ps.plotter(figure = fig2, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                          xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                          ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                          title = "CH$_3$CN Mean Spectrum")
plt.show()
########################################################
# Started Logging At: 2022-03-21 17:15:40
########################################################
freq_spw_2 = '135_spw47'
fmin = 135.25*u.GHz
fmax = 135.32*u.GHz
max_fn_2 = spectra+'source_ab_'+freq_spw_2+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.max.fits'
mean_fn_2 = spectra+'source_ab_'+freq_spw_2+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.mean.fits'
# Plot formaldehyde max spectrum
fig3 = plt.figure(3, figsize = (9, 5))
max_kspectrum = OneDSpectrum.from_hdu(fits.open(max_fn_2)).to(u.K)
max_kspectrum_ps = pyspeckit.Spectrum.from_hdu(max_kspectrum.hdu)
max_kspectrum_ps.xarr.convert_to_unit('GHz')
max_kspectrum_ps.plotter(figure = fig1, color = 'k', linewidth = 1, xmin = fmin, xmax = fmax,
                         xlabel = f"Frequency [{max_kspectrum_ps.xarr.unit.to_string('latex_inline')}]", 
                         ylabel = f"Brightness temperature [{max_kspectrum_ps.unit.to_string('latex_inline')}]",
                         title = "H$_2$CS Max Spectrum")
max_kspectrum_ps.plotter.savefig(f'../figures/h2cs_max_spec.pdf', facecolor='w', edgecolor='w', bbox_inches='tight')
plt.show()
########################################################
# # Started Logging At: 2022-03-21 17:15:43
########################################################
from pylab import imshow
from astropy.io import fits
from astropy.visualization import quantity_support
def generate_noise_map():
    '''Generate noise map from filename (hard-coded)'''
    hdu = fits.open('methyl_cyanide/template_noise.fits')
    noise_map = hdu[0].data*u.K # Generally, should NOT enforce units like this, but smarter solution is not implemented yet
    imshow(noise_map.value, origin='lower')
    return noise_map
noise_map = generate_noise_map()
import numpy as np
from lte_modeling_tools import nupper_of_kkms
from astropy import constants
def generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g):
    '''Get upper state column density maps from filenames (hard-coded)'''
    log_N_upper_gs = []
    log_N_upper_g_errs = []
    
    kmax = len(ch3cn_freqs) - 1 # Maximum rotational quantum number for the ladder

    for i in range(len(ch3cn_freqs)): # i is the index of the loop
        
        kk = kmax - i # Get rotational quantum number of this iteration

        if kk == 0 or kk == 1:
            # Import k = 0 and k = 1 cubes
            if kk == 0:
                # We need something to be called masked_cube to get channel_width and number_of_pixels
                masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
            elif kk == 1:
                # We need something to be called masked_cube to get channel_width and number_of_pixels
                masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_1_masked.fits', format='fits')
            masked_cube_0 = SpectralCube.read(f'methyl_cyanide/ch3cn_0_masked.fits', format='fits')
            print("k = 0:", masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]))
            masked_cube_1 = SpectralCube.read(f'methyl_cyanide/ch3cn_1_masked.fits', format='fits')
            print("k = 1:", masked_cube_1.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 1]))
            # Get frequency range we want: highest freq in k = 1, highest freq in k = 0
            highest_freq_1 = np.max(masked_cube_1.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 1]).spectral_axis)
            highest_freq_0 = np.max(masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]).spectral_axis)
            # Get extra part of k = 0 that's not in k = 1 cube
            extra_part_of_0 = masked_cube_0.with_spectral_unit(u.GHz, rest_value = ch3cn_freqs[kmax - 0]).spectral_slab(highest_freq_1, 
                                                                                                                        highest_freq_0)
            print("section of k = 0 that's not in k = 1:", extra_part_of_0)
            extra_part_of_0 = extra_part_of_0.with_spectral_unit(u.km/u.s, velocity_convention = 'radio')
            mom0_sum = extra_part_of_0.moment0() + masked_cube_1.moment0()
            # Split equally between k = 0 and k = 1 component
            mom0 = mom0_sum/2
        else: 
            # Import masked cube and get channel width
            masked_cube = SpectralCube.read(f'methyl_cyanide/ch3cn_{kk}_masked.fits', format='fits')
            # primary_beam = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.pb',
            #                                  format='casa_image')
            # masked_cube = masked_cube/primary_beam # Correct for effect of primary beam

            # Calculate moment 0 and moment 1 maps of cube
            mom0 = masked_cube.moment0()
            # mom1 = masked_cube.moment1()
    
        channel_width = np.diff(masked_cube.spectral_axis)[0]

        # Propagate error on integrated intensity
        number_of_pixels = masked_cube.mask.include().sum(axis=0)
        noise_map_int = noise_map*channel_width*np.sqrt(number_of_pixels) # Noise map WAS in Jy/beam... now K
            # Why were we not multiplying by the channel width???

        # Calculate what the shifted line frequency for the rung should be (per pixel) with mom1
        #shifted_line_freqs = (ch3cn_freqs[i]-((mom1/(constants.c.to(u.km/u.s)))*ch3cn_freqs[i])).to(u.GHz) # Maybe okay?

        # Calculate upper state column density from integrated line intensity (moment 0 map)
        N_upper = nupper_of_kkms(mom0, ch3cn_freqs[kmax - kk], ch3cn_A[kmax - kk])
        log_N_upper_g = np.log10(N_upper.value/ch3cn_g[kmax - kk]) # Shouldn't have to do .value?
        # Propagate error on upper state column density
        N_upper_err = nupper_of_kkms(noise_map_int, ch3cn_freqs[kmax - kk], ch3cn_A[kmax - kk])
        log_N_upper_g_err = N_upper_err/(N_upper*np.log(10.)) # There's no g here b/c it's a constant (divides out)
        # As of Feb. 23, 2022, I am once again convinced this is correct.

        # Append upper state column density maps and error maps into lists
        log_N_upper_gs.append(log_N_upper_g)
        log_N_upper_g_errs.append(log_N_upper_g_err)

    log_N_upper_gs = np.array(log_N_upper_gs)
    log_N_upper_g_errs = np.array(log_N_upper_g_errs)
    return log_N_upper_gs, log_N_upper_g_errs
log_N_upper_gs, log_N_upper_g_errs = generate_N_upper(ch3cn_freqs, ch3cn_A, ch3cn_g)
# Convert to natural log for fitting
ln_N_upper_gs = np.log(10**(log_N_upper_gs))
ln_N_upper_g_errs = np.log(10**(log_N_upper_g_errs))

# Replace all NaNs with 0s, for treatment as upper limits
ln_N_upper_gs = np.nan_to_num(ln_N_upper_gs)
# # Set the errors for NaN values to be huge so these are ignored in the fit
# ln_N_upper_g_errs = np.nan_to_num(ln_N_upper_g_errs, nan=1e10)
# Sum along k-component axis of NaN/non-NaN mask (which has True for "value," False for NaN)
# /where/ (not "if," as that invokes for loops) sum less than 3, make all values NaN (or 0)
ln_N_upper_gs.shape # (8, 512, 512)
ln_N_upper_gs_mask = ln_N_upper_gs != 0
# ln_N_upper_gs[ln_N_upper_gs_mask] # This does not have the shape I want, it's /just/ the valid values
ln_N_upper_gs_mask_sum = ln_N_upper_gs_mask.sum(axis=0)
ln_N_upper_gs_mask_sum_rep = np.repeat(ln_N_upper_gs_mask_sum[np.newaxis, :, :], 8, axis=0)

ln_N_upper_gs_test = ln_N_upper_gs.copy() # Probably don't need to duplicate this?
ln_N_upper_gs_test[ln_N_upper_gs_mask_sum_rep <= 3] = 0 # Use boolean mask to apply to upper state col densities

ln_N_upper_gs[ln_N_upper_gs_mask_sum_rep <= 3] = 0
# ln_N_upper_g_errs[ln_N_upper_gs_mask_sum_rep <= 3] = 1e10 # Set errors to be big so that they get ignored in the fit
# # This is my attempt at a brute-force solution to add upper limits. It's way too slow.
# flag = False
# for i in range(len(ch3cn_freqs)):
#     for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
#         while flag == False:
#             if ln_N_upper_gs[i,index[0],index[1]] == 0. and not np.isnan(noise_map[index[0],index[1]].value):
#                 ln_N_upper_g_errs[i,index[0],index[1]] = noise_map[index[0],index[1]].value
#                 flag = True
from scipy.optimize import curve_fit
import tqdm
print(tqdm.__version__) # As of Feb. 28, 2022: 4.62.3
from tqdm.notebook import tqdm
def linear(x, m, b):
    return m*x + b

def rotational_diagram(pixel_x, pixel_y, plot = False, save = False, print_results = False, plot_uncertainties = False):
    # Don't fit missing data
    col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pixel_y,pixel_x]) # ignores NaNs (shouldn't need this...)
    col_density_zero_mask = ln_N_upper_gs[:,pixel_y,pixel_x] != 0
    fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
#     fit_energies = ch3cn_E_U[::-1]
    fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x][col_density_zero_mask]
#     fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x]
    fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x][col_density_zero_mask]
#     fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x]
    
    # Do the simple linear fit
    # CONVERT TO K FROM ERGS TO DO THE FIT BECAUSE ERGS HAD SMALL NUMBERS AND THE FIT DIDN'T WORK
    fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()
    popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, sigma = fit_col_densities_errs)
    slope, intercept = popt[0], popt[1]
    
    # Extract the values
    temp = (-1./slope)*u.K
    # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BECAUSE FIT ALREADY DIVIDED BY K_B
    total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))

    if plot == True:
        plt.errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)), yerr = np.log10(np.exp(fit_col_densities_errs)), 
                     fmt = 'o') # NOW X AXIS IS IN KELVIN
        plt.plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)), 
                 label = f"$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {total_col_density:.2f}")
### Why were we exp'ing the y axis???
    #     plt.yscale('log')
### Why were we log scaling the y axis???

        ### PLOT UNCERTAINTIES:
        if plot_uncertainties == True:
            rng = np.random.default_rng()
            rv = rng.multivariate_normal(popt, pcov, 10)
            for param in rv:
                slope_i, intercept_i = param[0], param[1]
                temp_i = (-1./slope_i)*u.K
                total_col_density_i = np.log10(np.exp(intercept_i)*ch3cn_partfunc(temp_i))
                plt.plot(fit_energies_converted, np.log10(np.exp(slope_i*fit_energies_converted.value+intercept_i)),
                              alpha = 0.1, color = 'tab:orange') 
        plt.title(f"Rotational diagram for ({pixel_x}, {pixel_y})")
        plt.xlabel(f"Upper state energy [K]")
        plt.ylabel(f"Column density [$\log_{{10}}(N_{{upper}}\ [\mathrm{{cm}}^{{-2}}])$]")
        plt.legend()
        if save == True:
            plt.savefig(f"../figures/rot_diagram_x{pixel_x}y{pixel_y}.pdf", dpi = 200, facecolor='w', edgecolor='w', bbox_inches='tight')
        plt.show()

    if print_results == True:
        # Print extracted values
        print(f"Temp: {temp:.5f}")
        print(f"log10(Total column density): {total_col_density:.5f}")

    return temp.value, total_col_density
### Export a temperature map, for diagnostic purposes
### temp_map must already exist (for a single pixel)
# data = temp_map

# header = ch3cncube.header.copy() # or copy.copy(old_header)
# wcs = ch3cncube.wcs
# header.update(wcs.to_header())
# # header['BUNIT'] = temp_map.unit.to_string('fits')
# # header['RESTFRQ'] = ch3cn_freqs[comp].to(u.Hz).value

# fits.PrimaryHDU(data=data, header=header).writeto(f'temperature_map_februrary.fits', overwrite = True)
### Export a column density map, for diagnostic purposes
### temp_map must already exist (for a single pixel)
# data = total_col_density_map

# header = ch3cncube.header.copy() # or copy.copy(old_header)
# wcs = ch3cncube.wcs
# header.update(wcs.to_header())
# # header['BUNIT'] = temp_map.unit.to_string('fits')
# # header['RESTFRQ'] = ch3cn_freqs[comp].to(u.Hz).value

# fits.PrimaryHDU(data=data, header=header).writeto(f'col_density_map_februrary.fits', overwrite = True)
### Generate temperature map and total column density map
temp_map = np.zeros(ln_N_upper_gs.shape[1:])
total_col_density_map = np.zeros(ln_N_upper_gs.shape[1:])
counter = 0
for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
    col_density = ln_N_upper_gs[:,index[0],index[1]]
    col_density_errs = ln_N_upper_g_errs[:,index[0],index[1]]
    # Create column density mask, accounting for NaNs and zeros
    col_density_mask = (col_density != 0) & np.isfinite(col_density)
    if not col_density_mask.sum() > 1: # Are there unmasked values? If no, skip index
        continue
    temp_map[index], total_col_density_map[index] = rotational_diagram(index[1], index[0]) 
    # Should be x, y order for function
    counter += 1
    
    if counter % 30 == 0:
        print(f"x, y: {index[1]}, {index[0]}")
## Temperature map
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
plt.title('Temperature [K]')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar()
plt.savefig("../figures/temp_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
## Column density map
# mako_r = sns.color_palette("mako_r", as_cmap=True)
fig = plt.figure(dpi = 200)
# Set up normalization
from astropy.visualization import imshow_norm, ManualInterval, SqrtStretch, SinhStretch
cm = pl.matplotlib.cm.viridis.copy()
cm.set_under('w') # Make sure the "zero" color is white
# pl.imshow(total_col_density_map, cmap = cm, norm = norm, origin='lower')

plt.subplot(111,  projection = mom0.wcs)
im, norm = imshow_norm(total_col_density_map, origin='lower', cmap = cm, 
                       interval=ManualInterval(vmin = 5.0, vmax = np.max(total_col_density_map)), 
                       stretch=SinhStretch())
plt.tick_params(direction = 'in')
plt.title('Total Column Density [$\log_{{10}}(N_{{tot}}\ [\mathrm{{cm}}^{{-2}}])$]')
plt.xlabel('Right ascension'); plt.ylabel('Declination')
plt.colorbar(mappable=im)
plt.savefig("../figures/coldens_map.pdf", dpi = 300, facecolor='w', edgecolor='w', bbox_inches='tight')
# Get mom1/mom2 estimates
fn = '/blue/adamginsburg/abulatek/brick/first_results/temperature_map/methyl_cyanide/template_cube_masked.fits'
template_cube_masked = SpectralCube.read(fn, format='fits').with_spectral_unit(u.km/u.s)
h2cs_mom1 = template_cube_masked.moment1()
h2cs_linewidth_sigma = template_cube_masked.linewidth_sigma()
# Import statements

# Column density map
from astropy.visualization import imshow_norm, ManualInterval, SqrtStretch, SinhStretch
# Synthetic/real spectra
from astropy import units as u
import numpy as np
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule
plot_uncertainties = True

flag = False
for index in tqdm(np.ndindex(ln_N_upper_gs.shape[1:])):
    col_density = ln_N_upper_gs[:,index[0],index[1]]
    col_density_errs = ln_N_upper_g_errs[:,index[0],index[1]]
    # Create column density mask, accounting for NaNs and zeros
    col_density_mask = (col_density != 0) & np.isfinite(col_density)
    if not col_density_mask.sum() > 1: # Are there unmasked values? If no, skip index
        continue
    # For all valid pixels...
    
    ##################################################################################################  
    # Add plots: rotational diagram, synth/real spectra, and temp/col maps with X where the pixel is #
    ##################################################################################################
    
    if not flag:
        plt.clf()
        fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (10, 10), num = 0)
        
        ###### Temperature map
        
        cm = pl.matplotlib.cm.plasma.copy()
        cm.set_under('w') # Make sure the "zero" color is white

        ax[0][0].imshow(temp_map, vmax=200, vmin=0.001, cmap = cm, origin='lower')
        ax[0][0].scatter(index[1], index[0], marker = 'x', s = 100, linewidths = 2, color = 'r')
        ax[0][0].tick_params(direction = 'in')
        ax[0][0].set_title('Temperature (K)')
        # ax[0][0].colorbar()
        
        ###### Column density map

        cm = pl.matplotlib.cm.viridis.copy()
        cm.set_under('w') # Make sure the "zero" color is white

        im, norm = imshow_norm(total_col_density_map, origin='lower', cmap = cm, 
                               interval=ManualInterval(vmin = 10, 
                                                       vmax = np.max(total_col_density_map)), 
                               stretch=SinhStretch(), ax = ax[0][1])
        ax[0][1].scatter(index[1], index[0], marker = 'x', s = 100, linewidths = 2, color = 'r')
        ax[0][1].tick_params(direction = 'in')
        ax[0][1].set_title('Total Column Density ($\log_{10}$)')
        # ax[0][1].colorbar(mappable=im)
        
        ###### Rotational diagram
        
        ### Could use the built-in function, but don't know how to throw that plot to a subplot
        
        pixel_x, pixel_y = index[1], index[0]
        col_density_NaN_mask = np.isfinite(ln_N_upper_gs[:,pixel_y,pixel_x]) # ignores NaNs
        col_density_zero_mask = ln_N_upper_gs[:,pixel_y,pixel_x] != 0
        fit_energies = (ch3cn_E_U[::-1])[col_density_zero_mask]
        fit_col_densities = ln_N_upper_gs[:,pixel_y,pixel_x][col_density_zero_mask]
        fit_col_densities_errs = ln_N_upper_g_errs[:,pixel_y,pixel_x][col_density_zero_mask]

        # Do the simple linear fit
        # CONVERT TO K FROM ERGS TO DO THE FIT BECAUSE ERGS HAD SMALL NUMBERS AND THE FIT DIDN'T WORK
        fit_energies_converted = (fit_energies*u.erg/constants.k_B).decompose()
        popt, pcov = curve_fit(linear, fit_energies_converted, fit_col_densities, sigma = fit_col_densities_errs)
        slope, intercept = popt[0], popt[1]

        # Extract the values
        temp = (-1./slope)*u.K
        # temp = (temp*u.erg/constants.k_B).decompose() # DON'T NEED THIS BECAUSE FIT ALREADY DIVIDED BY K_B
        total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))
        
        ax[1][0].errorbar(fit_energies_converted, np.log10(np.exp(fit_col_densities)), yerr = np.log10(np.exp(fit_col_densities_errs)), 
                     fmt = 'o', color = 'tab:blue') 
        ax[1][0].plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)), 
                 label = f"$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {total_col_density:.2f}",
                 color = 'tab:orange') 
        ### PLOT UNCERTAINTIES:
        if plot_uncertainties == True:
            rng = np.random.default_rng()
            rv = rng.multivariate_normal(popt, pcov, 10)
            for param in rv:
                slope, intercept = param[0], param[1]
                temp = (-1./slope)*u.K
                total_col_density = np.log10(np.exp(intercept)*ch3cn_partfunc(temp))
                ax[1][0].plot(fit_energies_converted, np.log10(np.exp(slope*fit_energies_converted.value+intercept)),
                              alpha = 0.1, color = 'tab:orange') 
        ### Why were we exp'ing the y axis???
        # plt.yscale('log')
        ### Why were we log-scaling the y axis???
        ax[1][0].set_title(f"Rotational diagram for ({pixel_x}, {pixel_y})")
        ax[1][0].set_xlabel(f"Upper state energy [K]")
        ax[1][0].set_ylabel(f"Column density [$\log_{{10}}(N_{{upper}}\ [\mathrm{{cm}}^{{-2}}])$]")
        ax[1][0].legend()
        
        ###### Synthetic/real spectra
        
        fmin = 147.1*u.GHz
        fmax = 147.2*u.GHz
        
        sp_axis = np.linspace(fmin, fmax, 1000)
        
        fillingfactor = 1
        offset = 0
        species = 'CH3CN'

        freqs, aij, deg, EU, partfunc = ch3cn_freqs, ch3cn_A_orig, ch3cn_g, ch3cn_E_U, ch3cn_partfunc
        
#         # Get mom1 and mom2 values for a better estimation of the center velocity and velocity dispersion
#         masked_cube_2 = SpectralCube.read(f'methyl_cyanide/ch3cn_2_masked.fits', 
#                                           format='fits').with_spectral_unit(u.km/u.s, velocity_convention = 'radio')
#         mom1 = masked_cube_2.moment1()
#         linewidth_sigma = masked_cube_2.linewidth_sigma() # Yields nan...
        
        v_cen = h2cs_mom1[index[0], index[1]] # 10*u.km/u.s
        v_disp = h2cs_linewidth_sigma[index[0], index[1]] # 1.5*u.km/u.s
        # print(v_cen, v_disp)
        # temp = temperature
        N_tot = total_col_density # - np.log10(2) # This halves N_tot... do not run it unless you want to do that
        
        mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
        mod_sp = pyspeckit.Spectrum(xarr = sp_axis.to(u.GHz), data = mod, unit = u.K)
        
        # Get data spectrum
        data_sp = ch3cncube[:, pixel_y, pixel_x]
        data_sp_K = data_sp.value * ch3cncube.jtok_factors()
        data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

        mod_sp.plotter(label = f'$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {N_tot:.2f}',
                       color = 'tab:red', linewidth = 1, axis = ax[1][1])
        data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, #ymin = min([min(mod),min(data_sp)]) - 0.3, 
                            # ymax = max([max(mod),max(data_sp)]) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")

        plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
        plt.legend()
        if flag == False:
            plt.savefig(f"../diagnostic_plots/x{pixel_x}y{pixel_y}.png", overwrite = True)
        flag = True
        
print(counter)
from spectral_cube import SpectralCube
fig = plt.figure(dpi = 200) # I DEFINITELY shouldn't have to run this again
pixel_x, pixel_y = 339, 279 # 246, 274 # 173, 219 # 178, 226 # 227, 270
temperature, total_col_density = rotational_diagram(pixel_x, pixel_y, plot = True, save = True, plot_uncertainties = True)
fmin = 147.1*u.GHz
fmax = 147.2*u.GHz
sp_axis = np.linspace(fmin, fmax, 1000)
fillingfactor = 1
offset = 0
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
v_cen = h2cs_mom1[pixel_y, pixel_x] # 10*u.km/u.s
v_disp = h2cs_linewidth_sigma[pixel_y, pixel_x] # 1.5*u.km/u.s
# temp = 38.32444*u.K
# N_tot = (10**(13.71651))*u.cm**-2
temp = temperature*u.K
N_tot = 10**(total_col_density) / (u.cm**2) # - np.log10(2) # This halves N_tot, so don't run it if you don't want to half N_tot
mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

# mod_sp.plotter()
data_sp = ch3cncube[:, pixel_y, pixel_x]
data_sp_K = data_sp.value * ch3cncube.jtok_factors()
data_sp_K_pyspeckit = pyspeckit.Spectrum(xarr = ch3cncube.spectral_axis.to(u.GHz),
                                           data = data_sp_K, unit = u.K)

# data_sp_K_pyspeckit.plotter()
fig = plt.figure(dpi = 150)
mod_sp.plotter(label = f'$T$ = {temp:.1f}, $\log_{{10}}(N_{{tot}})$ = {np.log10(N_tot.value):.2f}',
               color = 'tab:red', linewidth = 1) # $v_{{c}}$ = {v_cen:.2f}, $v_{{disp}}$ = {v_disp:.2f}
data_sp_K_pyspeckit.plotter(axis=mod_sp.plotter.axis, clear = False, color = 'k', linewidth = 1,
                            xmin = fmin, xmax = fmax, #ymin = min([min(mod),min(data_sp)]) - 0.3, ymax = max([max(mod),max(data_sp)]) + 0.2, 
                            label = 'Data spectrum', 
                            xlabel = f"Frequency [{mod_sp.xarr.unit.to_string('latex_inline')}]", 
                            ylabel = f"Brightness temperature [{mod_sp.unit.to_string('latex_inline')}]")
                            
for i in range(len(ch3cn_freqs)): # i is the index of the loop
    if ch3cn_freqs[i].to(u.GHz).value > (plt.xlim()[0] + 0.01) and ch3cn_freqs[i].to(u.GHz).value < (plt.xlim()[1] - 0.01):
        kmax = len(ch3cn_freqs) - 1 # Maximum rotational quantum number for the ladder
        kk = kmax - i # Get rotational quantum number of this iteration
        shifted_freq = (ch3cn_freqs[kmax - kk]-(((v_cen)/(constants.c.to(u.km/u.s)))*ch3cn_freqs[kmax - kk])).to(u.GHz) # Maybe okay?

        plt.text(shifted_freq.value, plt.ylim()[1]-0.075, str(kk))
                            
plt.title(f"Synthetic and data spectra for ({pixel_x}, {pixel_y})")
plt.legend(loc = "lower left")
data_sp_K_pyspeckit.plotter.savefig(f'../figures/ch3cn_spec_datasynth_x{pixel_x}y{pixel_y}.pdf')
                            
print(f"Central velocity: {v_cen:.2f}")
print(f"Velocity dispersion: {v_disp:.2f}")
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these... but I don't want other lines creeping in
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
cube = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image')
from astropy import units as u
from spectral_cube import SpectralCube
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
freq_spw = '146_spw51'
fmin = 147.035*u.GHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these... but I don't want other lines creeping in
fmax = 147.175*u.GHz # ch3cncube.spectral_axis.max()
cube = SpectralCube.read(results+'source_ab_'+freq_spw+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image')
import regions
scube = cube.subcube_from_regionss([regions.Regions.read('../centralcoreellipse.reg')])
scube
scube = cube.subcube_from_regions([regions.Regions.read('../centralcoreellipse.reg')])
scube
scube = cube.subcube_from_regions(regions.Regions.read('../centralcoreellipse.reg'))
scube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 17, 17) and unit=Jy / beam and chunk size (128, 10, 17):
#[Out]#  n_x:     17  type_x: RA---SIN  unit_x: deg    range:   266.543777 deg:  266.544790 deg
#[Out]#  n_y:     17  type_y: DEC--SIN  unit_y: deg    range:   -28.705319 deg:  -28.704431 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
scube.spectral_slab(147.15*u.GHz, 147.155*u.GHz).moment0() # slab of comps 0 and 1
#[Out]# <Projection [[            nan,             nan,             nan,
#[Out]#                           nan,             nan,  66203.29539637,
#[Out]#                65123.08329578,  60836.23825048,  54269.81997343,
#[Out]#                46703.27878657,  39268.16061797,  32642.48414002,
#[Out]#                27047.23908462,             nan,             nan,
#[Out]#                           nan,             nan],
#[Out]#              [            nan,             nan,             nan,
#[Out]#                66315.95937074,  74580.57621738,  79522.50772397,
#[Out]#                80161.53545277,  76427.01962597,  69262.0538477 ,
#[Out]#                60227.87786893,  50875.74364795,  42257.4180688 ,
#[Out]#                34794.00742633,  28484.61829367,             nan,
#[Out]#                           nan,             nan],
#[Out]#              [            nan,             nan,  63591.49927719,
#[Out]#                78460.37485708,  91817.62358584, 101624.75348438,
#[Out]#               106130.01229876, 104519.31540349,  97269.94811053,
#[Out]#                86003.29464323,  72904.93496563,  59993.729839  ,
#[Out]#                48582.80920968,  39165.15314316,  31714.50338061,
#[Out]#                           nan,             nan],
#[Out]#              [            nan,             nan,  72496.39558109,
#[Out]#                93061.37736344, 113187.84944827, 129958.04666809,
#[Out]#               140492.17045082, 142797.05071019, 136485.54068083,
#[Out]#               123011.83004483, 105232.17214595,  86431.18067915,
#[Out]#                69270.8646063 ,  55184.03036972,  44449.78347262,
#[Out]#                36735.90345771,             nan],
#[Out]#              [            nan,  58880.32533312,  81071.66232656,
#[Out]#               107117.95097925, 134372.29925269, 159028.37393978,
#[Out]#               176894.47227235, 184541.94606031, 180472.7848934 ,
#[Out]#               165765.48892606, 143778.45948371, 118961.2603994 ,
#[Out]#                95356.61283333,  75545.2274051 ,  60437.79455053,
#[Out]#                49741.05623388,  42586.61952701],
#[Out]#              [            nan,  63580.73600411,  88055.56086511,
#[Out]#               118105.49605847, 151126.0677923 , 182756.60610046,
#[Out]#               207704.55157616, 221182.15587213, 220514.59919395,
#[Out]#               206178.89899964, 181662.72730795, 152152.93716059,
#[Out]#               122753.92532292,  97133.66523662,  77061.78761255,
#[Out]#                62667.42216827,  52955.21023621],
#[Out]#              [ 47578.90775117,  66279.65310616,  91625.63504373,
#[Out]#               123436.89747477, 159495.27115613, 195456.98093151,
#[Out]#               225561.99763419, 244154.10756814, 247554.13676627,
#[Out]#               235422.88750785, 210854.90903767, 179140.28095319,
#[Out]#               145918.96016971, 115677.27495616,  91069.48424483,
#[Out]#                72908.04348148,  60446.98526157],
#[Out]#              [ 46389.01905522,  64651.57639894,  89238.88014676,
#[Out]#               120269.28958069, 156013.32593718, 192655.81046135,
#[Out]#               224771.50189404, 246633.44853309, 254021.84340735,
#[Out]#               245730.58782363, 223959.63968343, 193394.95100027,
#[Out]#               159559.25232973, 127306.33462234,  99956.35997915,
#[Out]#                79009.51825042,  64192.11701904],
#[Out]#              [ 41478.51639243,  57755.55715652,  79707.47999268,
#[Out]#               107511.23155512, 139808.1923686 , 173519.42307212,
#[Out]#               204141.79813313, 226672.54112766, 237028.22686889,
#[Out]#               233403.52867007, 216875.63242761, 190945.43414277,
#[Out]#               160331.68861828, 129658.47405452, 102489.20125093,
#[Out]#                80788.34467468,  64754.93915545],
#[Out]#              [ 34721.81020726,  47533.0013205 ,  65247.96554966,
#[Out]#                88050.22743318, 114803.81292566, 143068.30447377,
#[Out]#               169398.86213744, 189962.07899403, 201439.45804829,
#[Out]#               201963.73422434, 191681.57690004, 172680.84731633,
#[Out]#               148357.85745846, 122521.54403779,  98503.53116806,
#[Out]#                78434.79548474,  62855.00567215],
#[Out]#              [            nan,  37667.69788186,  50635.35596199,
#[Out]#                67902.40365762,  88537.20028212, 110548.95255528,
#[Out]#               131323.4234306 , 148158.56864499, 158785.76554428,
#[Out]#               161791.2726702 , 156871.73471079, 144902.21195258,
#[Out]#               127819.50600012, 108284.42683905,  89069.03730734,
#[Out]#                72257.21312432,  58604.49895205],
#[Out]#              [            nan,  30892.6685408 ,  39940.44589651,
#[Out]#                52398.8587742 ,  67637.96789184,  84054.77385887,
#[Out]#                99590.31933382, 112344.9002601 , 120963.53251842,
#[Out]#               124673.68341097, 123147.79981609, 116486.17690286,
#[Out]#               105425.64388964,  91523.42487788,  76917.73310794,
#[Out]#                63563.86120377,  52375.16938678],
#[Out]#              [            nan,  27379.20196268,  34018.47660582,
#[Out]#                43012.7611014 ,  54089.95438166,  66116.24330971,
#[Out]#                77513.474907  ,  86912.36829154,  93541.87010168,
#[Out]#                97098.7976071 ,  97347.03093605,  93969.00601005,
#[Out]#                86921.0583112 ,  76958.66269582,  65701.22745257,
#[Out]#                54961.79519075,             nan],
#[Out]#              [            nan,             nan,  30735.6263018 ,
#[Out]#                37504.31275618,  45518.39188839,  54145.77612749,
#[Out]#                62365.56824832,  69304.87478477,  74578.53138649,
#[Out]#                78097.52208332,  79589.67199416,  78429.1762358 ,
#[Out]#                74080.41769952,  66787.17999579,  57790.13726368,
#[Out]#                48741.48955309,             nan],
#[Out]#              [            nan,             nan,             nan,
#[Out]#                32455.87871244,  38328.37107672,  44341.04690998,
#[Out]#                49984.15190533,  54936.27266348,  59236.73051157,
#[Out]#                62974.16430955,  65772.88566668,  66682.98173043,
#[Out]#                64739.96127513,  59763.20269797,  52649.93864849,
#[Out]#                           nan,             nan],
#[Out]#              [            nan,             nan,             nan,
#[Out]#                           nan,  30306.65458187,  34683.75720735,
#[Out]#                38435.01419173,  41642.48634021,  44838.18712512,
#[Out]#                48486.43708945,  52310.84560794,  55149.86841491,
#[Out]#                55604.00379918,  52982.64137403,             nan,
#[Out]#                           nan,             nan],
#[Out]#              [            nan,             nan,             nan,
#[Out]#                           nan,             nan,             nan,
#[Out]#                27955.87076251,  29931.00104781,  31840.6767708 ,
#[Out]#                34676.26937567,  38618.36955924,  42635.13639702,
#[Out]#                           nan,             nan,             nan,
#[Out]#                           nan,             nan]] Hz Jy / beam>
scube.spectral_slab(fmin,fmax).mean(axis=(1,2)).quicklook()
get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
scube.spectral_slab(fmin,fmax).mean(axis=(1,2)).quicklook()
scube.spectral_slab(147.15*u.GHz, 147.155*u.GHz).mean(axis=(1,2)).quicklook()
scube.spectral_slab(147.15*u.GHz, 147.157*u.GHz).mean(axis=(1,2)).quicklook()
scube.spectral_slab(147.15*u.GHz, 147.156*u.GHz).mean(axis=(1,2)).quicklook()
m0 = scube.spectral_slab(147.15*u.GHz, 147.156*u.GHz).moment0() # slab of comps 0 and 1
avgspec = (scube * m0).sum(axis=(1,2)) / m0.sum()
m0.sum()
#[Out]# <Projection nan Hz Jy / beam>
m0.value
#[Out]# array([[            nan,             nan,             nan,
#[Out]#                     nan,             nan,  71781.79540031,
#[Out]#          71785.80582005,  68149.66766443,  61866.25689255,
#[Out]#          54174.58924537,  46092.73401947,  38282.0149447 ,
#[Out]#          31156.35210419,             nan,             nan,
#[Out]#                     nan,             nan],
#[Out]#        [            nan,             nan,             nan,
#[Out]#          67931.0519248 ,  78735.44474395,  85700.49982367,
#[Out]#          87716.59370796,  84736.43948361,  77753.66920603,
#[Out]#          68337.45205878,  58065.13120553,  48130.64759996,
#[Out]#          39200.10433493,  31482.59037807,             nan,
#[Out]#                     nan,             nan],
#[Out]#        [            nan,             nan,  62217.39889958,
#[Out]#          80622.86816199,  97334.20782498, 109883.5248387 ,
#[Out]#         116272.01378079, 115597.0108531 , 108348.5239766 ,
#[Out]#          96255.89207595,  81737.25784979,  67137.9739973 ,
#[Out]#          54058.63075808,  43103.81473997,  34194.29137521,
#[Out]#                     nan,             nan],
#[Out]#        [            nan,             nan,  71967.39216896,
#[Out]#          96821.5037421 , 121188.44231482, 141608.90982962,
#[Out]#         154780.82202646, 158448.82252641, 152159.46115568,
#[Out]#         137531.13254164, 117793.91276359,  96690.85191248,
#[Out]#          77241.65814984,  61023.87251701,  48310.30481191,
#[Out]#          38791.17115617,             nan],
#[Out]#        [            nan,  55981.58609714,  82254.68328529,
#[Out]#         113095.92460214, 145327.24482072, 174532.59368788,
#[Out]#         195938.78247066, 205670.84686155, 202018.15724036,
#[Out]#         186144.53684849, 161776.07911968, 133891.89610228,
#[Out]#         107044.93985246,  84169.35317106,  66358.11226453,
#[Out]#          53411.20934297,  44530.55834286],
#[Out]#        [            nan,  62826.44352525,  91430.79503919,
#[Out]#         126466.27145842, 164877.64441051, 201698.1458679 ,
#[Out]#         230990.08978735, 247386.41106467, 247798.10034826,
#[Out]#         232566.47272409, 205411.45916543, 172097.50229031,
#[Out]#         138474.83939041, 108871.08547793,  85445.82526606,
#[Out]#          68426.34990488,  56723.82707991],
#[Out]#        [ 46191.48948299,  67858.45195841,  97203.80875697,
#[Out]#         133908.56171276, 175371.22083162, 216729.66610958,
#[Out]#         251647.46845728, 273841.76230392, 279008.85190618,
#[Out]#         266389.23709459, 239130.58405763, 203150.29367503,
#[Out]#         165070.30645301, 130286.35355542, 101948.48980064,
#[Out]#          80905.73793156,  66191.5931564 ],
#[Out]#        [ 47147.68803607,  68129.55380353,  96394.90010255,
#[Out]#         131915.95507482, 172649.17688097, 214432.76053151,
#[Out]#         251432.17298139, 277306.56123048, 287031.357638  ,
#[Out]#         278731.27242313, 254518.55472644, 219741.36434585,
#[Out]#         181025.48633341, 144201.6442049 , 113052.05862529,
#[Out]#          89052.62226994,  71716.32841024],
#[Out]#        [ 43757.91277918,  62252.94076849,  87231.81962636,
#[Out]#         118719.73683523, 155180.43059561, 193421.69683346,
#[Out]#         228708.61150673, 255429.73482569, 268558.06936388,
#[Out]#         265468.86719381, 247086.04262011, 217563.42504146,
#[Out]#         182672.83556605, 147889.18121878, 117116.85409077,
#[Out]#          92294.07714073,  73540.7616348 ],
#[Out]#        [ 37860.03556639,  52151.04314284,  71906.31884811,
#[Out]#          97245.24324276, 127081.75111665, 159095.47749943,
#[Out]#         189711.51046957, 214443.60293293, 228943.36254917,
#[Out]#         230492.04874837, 219089.12218451, 197420.06496824,
#[Out]#         169762.85977092, 140561.48823967, 113368.58766872,
#[Out]#          90335.82880048,  72068.95702957],
#[Out]#        [            nan,  41903.86355696,  55904.10865314,
#[Out]#          74578.12835494,  97278.85505098, 122279.34674584,
#[Out]#         146835.19811599, 167535.55711657, 181088.13869664,
#[Out]#         185309.04073124, 179799.54916528, 165953.7970909 ,
#[Out]#         146401.07049129, 124262.07916761, 102470.06538106,
#[Out]#          83166.27062263,  67205.27094388],
#[Out]#        [            nan,  34729.83474168,  44221.67763317,
#[Out]#          57419.61606145,  74078.27151387,  92852.64294562,
#[Out]#         111501.86603171, 127442.56542616, 138441.60094062,
#[Out]#         143128.82992182, 141142.99694252, 133013.8935149 ,
#[Out]#         120023.82460669, 104105.33412664,  87552.53678722,
#[Out]#          72362.80284332,  59452.52872526],
#[Out]#        [            nan,  31023.1917016 ,  38104.23954619,
#[Out]#          47774.8895458 ,  60045.49500689,  73912.17675419,
#[Out]#          87593.28083278,  99197.45236269, 107350.93755069,
#[Out]#         111373.63355673, 111039.47367353, 106353.14966585,
#[Out]#          97688.30509355,  86120.8748067 ,  73417.52208977,
#[Out]#          61380.57914486,             nan],
#[Out]#        [            nan,             nan,  35039.14909399,
#[Out]#          42737.27982676,  51906.2082396 ,  61899.70356584,
#[Out]#          71539.8535701 ,  79670.81537677,  85613.20102941,
#[Out]#          89116.79239656,  89954.38746119,  87708.91484841,
#[Out]#          82102.18671922,  73573.63582138,  63437.43452615,
#[Out]#          53292.33331084,             nan],
#[Out]#        [            nan,             nan,             nan,
#[Out]#          37743.81292005,  44652.57222289,  51605.53749296,
#[Out]#          58003.61437944,  63408.88221346,  67772.66842635,
#[Out]#          71197.49981102,  73428.12426392,  73657.49366982,
#[Out]#          70987.968539  ,  65204.15343352,  57122.49173666,
#[Out]#                     nan,             nan],
#[Out]#        [            nan,             nan,             nan,
#[Out]#                     nan,  35360.43217595,  40410.49915119,
#[Out]#          44642.13541105,  47992.60110774,  50935.41352085,
#[Out]#          54079.45123162,  57424.97174822,  59988.881125  ,
#[Out]#          60277.47701224,  57294.5117796 ,             nan,
#[Out]#                     nan,             nan],
#[Out]#        [            nan,             nan,             nan,
#[Out]#                     nan,             nan,             nan,
#[Out]#          32051.52707585,  34203.45474022,  35731.91116077,
#[Out]#          37869.859631  ,  41277.01430192,  45216.89127477,
#[Out]#                     nan,             nan,             nan,
#[Out]#                     nan,             nan]])
np.nansum(m0.value)
#[Out]# 26140532.02082247
avgspec = (scube * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
m0.quicklook()
avgspec
#[Out]# <VaryingResolutionOneDSpectrum [ 0.00870287,-0.00484011,-0.01063703,...,
#[Out]#                                  0.00732376, 0.00634538, 0.00523109] Jy / beam>
avgspec.quicklook()
avgspec.spectral_slab(fmin, fmax).quicklook()
avgspec = (scube.spectral_slab(fmin, fmax) * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
avgspec.quicklook()
meanspec = scube.spectral_slab(fmin, fmax).mean(axis=(1,2))
avgspec.quicklook()
pl.plot(meanspec.spectral_axis, meanspec)
avgspec.quicklook()
pl.plot(meanspec.spectral_axis, meanspec.value)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3ccfe6d250>]
reg = regions.Regions.read('../centralcoreellipse.reg')
scube = cube.subcube_from_regions(reg)
scube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 17, 17) and unit=Jy / beam and chunk size (128, 10, 17):
#[Out]#  n_x:     17  type_x: RA---SIN  unit_x: deg    range:   266.543777 deg:  266.544790 deg
#[Out]#  n_y:     17  type_y: DEC--SIN  unit_y: deg    range:   -28.705319 deg:  -28.704431 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 145582599820.702 Hz:147457166014.871 Hz
tmap = fits.open('temperature_map_februrary.fits')
from astropy.io import fits
tmap = fits.open('temperature_map_februrary.fits')
pl.imshow(tmap[0].data)
#[Out]# <matplotlib.image.AxesImage at 0x2b3ccfef2f10>
tmap[0].max()
tmap[0].data.max()
#[Out]# 40793.157454977794
pl.imshow(tmap[0].data, vmin=50, vmax=200)
#[Out]# <matplotlib.image.AxesImage at 0x2b3cd0169a90>
import pyspeckit
sp = pyspeckit.Spectrum(data=avgspec, xarr=avgspec.spectral_axis)
sp.plotter()
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule
sp = pyspeckit.Spectrum(data=avgspec, xarr=avgspec.spectral_axis)
sp.plotter()
vcen = 0*u.km/u.s
vdisp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e13*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
v_cen = 0*u.km/u.s
vdisp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e13*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
v_cen = 0*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e13*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
pl.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d02f6f8e0>]
v_cen = 10*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e13*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
pl.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d02f140d0>]
v_cen = 10*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e17*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
pl.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d02fc2550>]
v_cen = 20*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e17*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
sp.plotter()
sp.plotter.ax.plot(sp.xarr, mod)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d02633340>]
avgspec.beams
#[Out]# <Beams [4.78156366e-11, 4.78153493e-11, 4.78150620e-11, 4.78147861e-11,
#[Out]#         4.78145026e-11, 4.78142067e-11, 4.78139291e-11, 4.78136483e-11,
#[Out]#         4.78133696e-11, 4.78130899e-11, 4.78128092e-11, 4.78125144e-11,
#[Out]#         4.78122175e-11, 4.78119378e-11, 4.78116495e-11, 4.78113509e-11,
#[Out]#         4.78110846e-11, 4.78108059e-11, 4.78105063e-11, 4.78102180e-11,
#[Out]#         4.78099325e-11, 4.78096394e-11, 4.78093419e-11, 4.78090519e-11,
#[Out]#         4.78087647e-11, 4.78084668e-11, 4.78081806e-11, 4.78078820e-11,
#[Out]#         4.78075862e-11, 4.78073038e-11, 4.78070251e-11, 4.78067330e-11,
#[Out]#         4.78064506e-11, 4.78061795e-11, 4.78058789e-11, 4.78055680e-11,
#[Out]#         4.78052807e-11, 4.78050086e-11, 4.78047251e-11, 4.78044170e-11,
#[Out]#         4.78041335e-11, 4.78038387e-11, 4.78035515e-11, 4.78032595e-11,
#[Out]#         4.78029695e-11, 4.78026871e-11, 4.78024047e-11, 4.78021222e-11,
#[Out]#         4.78018264e-11, 4.78015334e-11, 4.78012499e-11, 4.78009685e-11,
#[Out]#         4.78006765e-11, 4.78003817e-11, 4.78000688e-11, 4.77997987e-11,
#[Out]#         4.77994992e-11, 4.77992178e-11, 4.77989446e-11, 4.77986536e-11,
#[Out]#         4.77983482e-11, 4.77980648e-11, 4.77977920e-11, 4.77975106e-11,
#[Out]#         4.77972309e-11, 4.77969372e-11, 4.77966671e-11, 4.77963772e-11,
#[Out]#         4.77960900e-11, 4.77958161e-11, 4.77955108e-11, 4.77952225e-11,
#[Out]#         4.77949316e-11, 4.77946293e-11, 4.77943431e-11, 4.77940377e-11,
#[Out]#         4.77937516e-11, 4.77934681e-11, 4.77931895e-11, 4.77929050e-11,
#[Out]#         4.77926226e-11, 4.77923241e-11, 4.77920369e-11, 4.77917573e-11,
#[Out]#         4.77914776e-11, 4.77911819e-11, 4.77908916e-11, 4.77905797e-11,
#[Out]#         4.77903000e-11, 4.77900118e-11, 4.77897133e-11, 4.77894337e-11,
#[Out]#         4.77891513e-11, 4.77888593e-11, 4.77885635e-11, 4.77882736e-11,
#[Out]#         4.77879827e-11, 4.77877078e-11, 4.77874158e-11, 4.77871040e-11,
#[Out]#         4.77868243e-11, 4.77865286e-11, 4.77862424e-11, 4.77859429e-11,
#[Out]#         4.77856606e-11, 4.77853754e-11, 4.77850920e-11, 4.77848059e-11,
#[Out]#         4.77845273e-11, 4.77842476e-11, 4.77839567e-11, 4.77836572e-11,
#[Out]#         4.77833690e-11, 4.77830781e-11, 4.77827813e-11, 4.77824893e-11,
#[Out]#         4.77821936e-11, 4.77819044e-11, 4.77816162e-11, 4.77813290e-11,
#[Out]#         4.77810494e-11, 4.77807708e-11, 4.77804912e-11, 4.77802098e-11,
#[Out]#         4.77799227e-11, 4.77796602e-11, 4.77793692e-11, 4.77790811e-11,
#[Out]#         4.77787874e-11, 4.77785050e-11, 4.77782237e-11, 4.77779280e-11,
#[Out]#         4.77776377e-11, 4.77773667e-11, 4.77770795e-11, 4.77768058e-11,
#[Out]#         4.77765176e-11, 4.77762246e-11, 4.77759385e-11, 4.77756352e-11,
#[Out]#         4.77753416e-11, 4.77750459e-11, 4.77747577e-11, 4.77744753e-11,
#[Out]#         4.77741882e-11, 4.77738946e-11, 4.77736283e-11, 4.77733305e-11,
#[Out]#         4.77730461e-11, 4.77727600e-11, 4.77724804e-11, 4.77721789e-11,
#[Out]#         4.77718907e-11, 4.77716009e-11, 4.77713137e-11, 4.77710218e-11,
#[Out]#         4.77707206e-11, 4.77704297e-11, 4.77701426e-11, 4.77698603e-11,
#[Out]#         4.77695666e-11, 4.77692559e-11, 4.77689660e-11, 4.77686693e-11,
#[Out]#         4.77683832e-11, 4.77680933e-11, 4.77678130e-11, 4.77675335e-11,
#[Out]#         4.77672501e-11, 4.77669705e-11, 4.77666796e-11, 4.77663877e-11,
#[Out]#         4.77660931e-11, 4.77658049e-11, 4.77655130e-11, 4.77652392e-11,
#[Out]#         4.77649597e-11, 4.77646640e-11, 4.77643683e-11, 4.77640774e-11,
#[Out]#         4.77637930e-11, 4.77634813e-11, 4.77631856e-11, 4.77629070e-11,
#[Out]#         4.77626151e-11, 4.77623232e-11, 4.77620372e-11, 4.77617473e-11,
#[Out]#         4.77614554e-11, 4.77611930e-11, 4.77609059e-11, 4.77606198e-11,
#[Out]#         4.77603118e-11, 4.77600322e-11, 4.77597489e-11, 4.77594361e-11,
#[Out]#         4.77591613e-11, 4.77588893e-11, 4.77585937e-11, 4.77583103e-11,
#[Out]#         4.77580147e-11, 4.77577330e-11, 4.77574439e-11, 4.77571455e-11,
#[Out]#         4.77568680e-11, 4.77565847e-11, 4.77562955e-11, 4.77560095e-11,
#[Out]#         4.77557101e-11, 4.77554401e-11, 4.77551605e-11, 4.77548755e-11,
#[Out]#         4.77546083e-11, 4.77543212e-11, 4.77540331e-11, 4.77537364e-11,
#[Out]#         4.77534514e-11, 4.77531595e-11, 4.77528649e-11, 4.77525837e-11,
#[Out]#         4.77522928e-11, 4.77520095e-11, 4.77517337e-11, 4.77514467e-11,
#[Out]#         4.77511596e-11, 4.77508708e-11, 4.77505601e-11, 4.77502549e-11,
#[Out]#         4.77499555e-11, 4.77496732e-11, 4.77493786e-11, 4.77490963e-11,
#[Out]#         4.77488141e-11, 4.77485232e-11, 4.77482410e-11, 4.77479491e-11,
#[Out]#         4.77476668e-11, 4.77473722e-11, 4.77470975e-11, 4.77468056e-11,
#[Out]#         4.77465271e-11, 4.77462401e-11, 4.77459654e-11, 4.77456677e-11,
#[Out]#         4.77453940e-11, 4.77451097e-11, 4.77448349e-11, 4.77445602e-11,
#[Out]#         4.77442598e-11, 4.77439728e-11, 4.77436819e-11, 4.77433798e-11,
#[Out]#         4.77430965e-11, 4.77428198e-11, 4.77425289e-11, 4.77422477e-11,
#[Out]#         4.77419483e-11, 4.77416592e-11, 4.77413732e-11, 4.77410841e-11,
#[Out]#         4.77407998e-11, 4.77405117e-11, 4.77402069e-11, 4.77399284e-11,
#[Out]#         4.77396318e-11, 4.77393506e-11, 4.77390635e-11, 4.77387792e-11,
#[Out]#         4.77384960e-11, 4.77381966e-11, 4.77379229e-11, 4.77376342e-11,
#[Out]#         4.77373499e-11, 4.77370526e-11, 4.77367591e-11, 4.77364635e-11,
#[Out]#         4.77361545e-11, 4.77358760e-11, 4.77355805e-11, 4.77352838e-11,
#[Out]#         4.77350102e-11, 4.77347108e-11, 4.77344296e-11, 4.77341388e-11,
#[Out]#         4.77338566e-11, 4.77335819e-11, 4.77332911e-11, 4.77330127e-11] sr>
avgspec.beams.common_beam()
#[Out]# Beam: BMAJ=1.512071132659912 arcsec BMIN=1.1873631477355957 arcsec BPA=-81.83331298828125 deg
avgspec.beams.common_beam().jtok(avgspec.spectral_axis)
#[Out]# <Quantity [31.48603485, 31.48582572, 31.4856166 , 31.48540748, 31.48519836,
#[Out]#            31.48498924, 31.48478012, 31.48457101, 31.4843619 , 31.48415279,
#[Out]#            31.48394368, 31.48373457, 31.48352547, 31.48331637, 31.48310727,
#[Out]#            31.48289817, 31.48268908, 31.48247999, 31.4822709 , 31.48206181,
#[Out]#            31.48185272, 31.48164364, 31.48143455, 31.48122547, 31.4810164 ,
#[Out]#            31.48080732, 31.48059824, 31.48038917, 31.4801801 , 31.47997104,
#[Out]#            31.47976197, 31.47955291, 31.47934384, 31.47913478, 31.47892573,
#[Out]#            31.47871667, 31.47850762, 31.47829857, 31.47808952, 31.47788047,
#[Out]#            31.47767143, 31.47746238, 31.47725334, 31.4770443 , 31.47683527,
#[Out]#            31.47662623, 31.4764172 , 31.47620817, 31.47599914, 31.47579012,
#[Out]#            31.47558109, 31.47537207, 31.47516305, 31.47495403, 31.47474502,
#[Out]#            31.474536  , 31.47432699, 31.47411798, 31.47390897, 31.47369997,
#[Out]#            31.47349097, 31.47328196, 31.47307296, 31.47286397, 31.47265497,
#[Out]#            31.47244598, 31.47223699, 31.472028  , 31.47181901, 31.47161003,
#[Out]#            31.47140105, 31.47119207, 31.47098309, 31.47077411, 31.47056514,
#[Out]#            31.47035617, 31.4701472 , 31.46993823, 31.46972926, 31.4695203 ,
#[Out]#            31.46931134, 31.46910238, 31.46889342, 31.46868446, 31.46847551,
#[Out]#            31.46826656, 31.46805761, 31.46784866, 31.46763972, 31.46743078,
#[Out]#            31.46722183, 31.4670129 , 31.46680396, 31.46659502, 31.46638609,
#[Out]#            31.46617716, 31.46596823, 31.46575931, 31.46555038, 31.46534146,
#[Out]#            31.46513254, 31.46492362, 31.46471471, 31.46450579, 31.46429688,
#[Out]#            31.46408797, 31.46387906, 31.46367016, 31.46346126, 31.46325235,
#[Out]#            31.46304345, 31.46283456, 31.46262566, 31.46241677, 31.46220788,
#[Out]#            31.46199899, 31.4617901 , 31.46158122, 31.46137234, 31.46116346,
#[Out]#            31.46095458, 31.4607457 , 31.46053683, 31.46032795, 31.46011908,
#[Out]#            31.45991022, 31.45970135, 31.45949249, 31.45928362, 31.45907476,
#[Out]#            31.45886591, 31.45865705, 31.4584482 , 31.45823935, 31.4580305 ,
#[Out]#            31.45782165, 31.45761281, 31.45740396, 31.45719512, 31.45698628,
#[Out]#            31.45677745, 31.45656861, 31.45635978, 31.45615095, 31.45594212,
#[Out]#            31.45573329, 31.45552447, 31.45531565, 31.45510683, 31.45489801,
#[Out]#            31.45468919, 31.45448038, 31.45427157, 31.45406276, 31.45385395,
#[Out]#            31.45364514, 31.45343634, 31.45322754, 31.45301874, 31.45280994,
#[Out]#            31.45260115, 31.45239235, 31.45218356, 31.45197477, 31.45176599,
#[Out]#            31.4515572 , 31.45134842, 31.45113964, 31.45093086, 31.45072208,
#[Out]#            31.45051331, 31.45030454, 31.45009576, 31.449887  , 31.44967823,
#[Out]#            31.44946947, 31.4492607 , 31.44905194, 31.44884319, 31.44863443,
#[Out]#            31.44842568, 31.44821693, 31.44800818, 31.44779943, 31.44759068,
#[Out]#            31.44738194, 31.4471732 , 31.44696446, 31.44675572, 31.44654699,
#[Out]#            31.44633826, 31.44612952, 31.4459208 , 31.44571207, 31.44550334,
#[Out]#            31.44529462, 31.4450859 , 31.44487718, 31.44466847, 31.44445975,
#[Out]#            31.44425104, 31.44404233, 31.44383362, 31.44362492, 31.44341621,
#[Out]#            31.44320751, 31.44299881, 31.44279011, 31.44258142, 31.44237273,
#[Out]#            31.44216403, 31.44195534, 31.44174666, 31.44153797, 31.44132929,
#[Out]#            31.44112061, 31.44091193, 31.44070325, 31.44049458, 31.44028591,
#[Out]#            31.44007724, 31.43986857, 31.4396599 , 31.43945124, 31.43924257,
#[Out]#            31.43903391, 31.43882526, 31.4386166 , 31.43840795, 31.43819929,
#[Out]#            31.43799064, 31.437782  , 31.43757335, 31.43736471, 31.43715607,
#[Out]#            31.43694743, 31.43673879, 31.43653015, 31.43632152, 31.43611289,
#[Out]#            31.43590426, 31.43569563, 31.43548701, 31.43527839, 31.43506977,
#[Out]#            31.43486115, 31.43465253, 31.43444392, 31.4342353 , 31.43402669,
#[Out]#            31.43381808, 31.43360948, 31.43340087, 31.43319227, 31.43298367,
#[Out]#            31.43277507, 31.43256648, 31.43235789, 31.43214929, 31.4319407 ,
#[Out]#            31.43173212, 31.43152353, 31.43131495, 31.43110637, 31.43089779,
#[Out]#            31.43068921, 31.43048064, 31.43027206, 31.43006349, 31.42985492,
#[Out]#            31.42964636, 31.42943779, 31.42922923, 31.42902067, 31.42881211,
#[Out]#            31.42860355, 31.428395  , 31.42818645, 31.4279779 , 31.42776935,
#[Out]#            31.4275608 , 31.42735226, 31.42714372, 31.42693518, 31.42672664,
#[Out]#            31.4265181 , 31.42630957, 31.42610104] K>
jtok = avgspec.beams.common_beam().jtok(avgspec.spectral_axis).mean()
jtok
#[Out]# <Quantity 31.45605372 K>
sp = pyspeckit.Spectrum(data=avgspec*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 20*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e17*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0922ddc0>]
v_cen = 20*u.km/u.s
v_disp = 5*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d09299b20>]
v_cen = 20*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
v_cen = 20*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d092dafa0>]
ff = 1.2/20
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod/ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0937ee80>]
ff = 1.2/20 / 10
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod/ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d093e3d60>]
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d09452e50>]
ff = 1.2/20
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d094c8d30>]
ff = 1.2/20 / 5
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d09532c40>]
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
ff = 1.2/20 / 5
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d095d2c70>]
avgspec = (scube.spectral_slab(147.000*u.GHz, 147.195*u.GHz) * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
avgspec.quicklook()
pl.plot(meanspec.spectral_axis, meanspec.value)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d096344f0>]
sp = pyspeckit.Spectrum(data=avgspec*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec.quantity))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
ff = 1.2/20 / 5
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3cfc212f10>]
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
fitter = generate_fitter(mod, name="CH3CN")
sp.specfit(fitter=fitter)
sp.plotter(plot_fit=True)
sp.plotter()
get_ipython().run_line_magic('pinfo', 'sp.plotter')
get_ipython().run_line_magic('pinfo2', 'sp.plotter')
sp.plotter()
sp.specfit.plot_fit()
get_ipython().run_line_magic('pinfo', 'sp.specfit')
sp.specfit(fittype=fitter)
get_ipython().run_line_magic('pinfo', 'sp.specfit.register_fitter')
fitter
#[Out]# <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3cfc37d520>
sp.specfit.register_fitter('CH3CN', function=fitter)
get_ipython().run_line_magic('pinfo', 'fitter')
get_ipython().run_line_magic('pinfo', 'fitter.modelfunc')
sp.specfit.register_fitter('CH3CN', function=fitter, npars=4)
sp.specfit(fittype=fitter)
sp.specfit(fittype='CH3CN')
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=4)
sp.specfit(fittype='CH3CN', guesses=[v_cen, v_disp, temp, N_tot])
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=4)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
get_ipython().run_line_magic('debug', '')
fitter.modelfunc
#[Out]# array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 3.61758342e-14, 1.93540703e-12, 6.98012651e-11,
#[Out]#        1.96145344e-09, 4.30258826e-08, 7.36652525e-07, 9.84413438e-06,
#[Out]#        1.02677275e-04, 8.35896187e-04, 5.31138061e-03, 2.63404079e-02,
#[Out]#        1.01940525e-01, 3.07807651e-01, 7.24955302e-01, 1.33189302e+00,
#[Out]#        1.91005425e+00, 2.14028244e+00, 1.87476215e+00, 1.28304191e+00,
#[Out]#        6.85366883e-01, 2.85575531e-01, 9.28162276e-02, 2.35366198e-02,
#[Out]#        4.65778169e-03, 7.19409735e-04, 8.67263651e-05, 8.16032206e-06,
#[Out]#        5.99301574e-07, 3.43530668e-08, 1.53698236e-09, 5.36667698e-11,
#[Out]#        1.46511902e-12, 3.61757765e-14, 0.00000000e+00, 0.00000000e+00,
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
#[Out]#        9.04392474e-14, 4.37725936e-12, 1.68904322e-10, 5.08767719e-09,
#[Out]#        1.19626013e-07, 2.19567296e-06, 3.14590411e-05, 3.51850993e-04,
#[Out]#        3.07188352e-03, 2.09346826e-02, 1.11344084e-01, 4.61908953e-01,
#[Out]#        1.49268661e+00, 3.75099566e+00, 7.32666065e+00, 1.11594661e+01,
#[Out]#        1.33377158e+01, 1.25648527e+01, 9.31442518e+00, 5.40034875e+00,
#[Out]#        2.43631429e+00, 8.54231656e-01, 2.33108148e-01, 4.95894983e-02,
#[Out]#        8.23129544e-03, 1.06644127e-03, 1.07853261e-04, 8.51460171e-06,
#[Out]#        5.24724236e-07, 2.52425885e-08, 9.47928510e-10, 2.77828939e-11,
#[Out]#        6.33073722e-13, 1.80878197e-14, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        2.71316877e-13, 1.23177856e-11, 4.46732237e-10, 1.26490622e-08,
#[Out]#        2.79605316e-07, 4.82518452e-06, 6.50075414e-05, 6.83745058e-04,
#[Out]#        5.61436058e-03, 3.59876822e-02, 1.80026832e-01, 7.02257606e-01,
#[Out]#        2.13262396e+00, 5.03264483e+00, 9.23182589e+00, 1.32257547e+01,
#[Out]#        1.49030200e+01, 1.32552298e+01, 9.27357304e+00, 5.06728983e+00,
#[Out]#        2.15240233e+00, 7.10444702e-01, 1.82552425e-01, 3.65776485e-02,
#[Out]#        5.71968103e-03, 6.98191868e-04, 6.65354809e-05, 4.95008257e-06,
#[Out]#        2.87509969e-07, 1.30369199e-08, 4.61509317e-10, 1.27518735e-11,
#[Out]#        2.71316444e-13, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        2.17052939e-13, 1.06717690e-11, 4.06087920e-10, 1.20659349e-08,
#[Out]#        2.79916913e-07, 5.07009564e-06, 7.17004250e-05, 7.91671628e-04,
#[Out]#        6.82464388e-03, 4.59292618e-02, 2.41217730e-01, 9.87439226e-01,
#[Out]#        3.14213837e+00, 7.74549928e+00, 1.47824817e+01, 2.19938315e+01,
#[Out]#        2.58318646e+01, 2.41512186e+01, 1.79063291e+01, 1.04003163e+01,
#[Out]#        4.68530445e+00, 1.63364442e+00, 4.42138245e-01, 9.31850910e-02,
#[Out]#        1.53211547e-02, 1.96638793e-03, 1.97038050e-04, 1.54151420e-05,
#[Out]#        9.41593344e-07, 4.49053302e-08, 1.67206481e-09, 4.86016954e-11,
#[Out]#        1.10335068e-12, 1.80877151e-14, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        5.42631265e-14, 3.05682264e-12, 1.35476926e-10, 4.69055821e-09,
#[Out]#        1.26806552e-07, 2.67674112e-06, 4.41182851e-05, 5.67777077e-04,
#[Out]#        5.70530096e-03, 4.47591339e-02, 2.74017724e-01, 1.30661364e+00,
#[Out]#        4.82700837e+00, 1.36814076e+01, 2.94972065e+01, 4.87006446e+01,
#[Out]#        6.34143202e+01, 6.74520501e+01, 5.92848623e+01, 4.21023527e+01,
#[Out]#        2.32926996e+01, 9.81949087e+00, 3.15538419e+00, 7.80856587e-01,
#[Out]#        1.50066428e-01, 2.24838563e-02, 2.62943246e-03, 2.40091940e-04,
#[Out]#        1.71174031e-05, 9.52896076e-07, 4.14209721e-08, 1.48897789e-09,
#[Out]#        2.93671572e-09, 7.88674538e-08, 1.67509034e-06, 2.77812138e-05,
#[Out]#        3.59775940e-04, 3.63812329e-03, 2.87251272e-02, 1.77031478e-01,
#[Out]#        8.50560182e-01, 3.17466589e+00, 9.14474865e+00, 2.02059074e+01,
#[Out]#        3.43789382e+01, 4.59400234e+01, 4.94268927e+01, 4.31970376e+01,
#[Out]#        3.01982462e+01, 1.64619762e+01, 6.89118975e+00, 2.21555061e+00,
#[Out]#        5.50933077e-01, 1.06718613e-01, 1.77700659e-02, 1.67101829e-02,
#[Out]#        1.01509902e-01, 5.40632128e-01, 2.24338737e+00, 7.19495709e+00,
#[Out]#        1.76845147e+01, 3.32590130e+01, 4.86690362e+01, 5.72983288e+01,
#[Out]#        5.63250082e+01, 4.90071031e+01, 4.39471675e+01, 4.76690915e+01,
#[Out]#        5.63346572e+01, 6.06864869e+01, 5.53505801e+01, 4.11952097e+01,
#[Out]#        2.41160954e+01, 1.08297282e+01, 3.71491613e+00, 9.81103683e-01,
#[Out]#        2.01105257e-01, 3.21292172e-02, 4.00658221e-03, 3.90121811e-04,
#[Out]#        2.96624459e-05, 1.76115723e-06, 8.16533718e-08, 2.95620672e-09,
#[Out]#        8.35829523e-11, 1.84493848e-12, 3.61752625e-14, 0.00000000e+00,
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
fitter = generate_fitter(mod, name="CH3CN")
fitter.modelfunc
#[Out]# array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 3.61758342e-14, 1.93540703e-12, 6.98012651e-11,
#[Out]#        1.96145344e-09, 4.30258826e-08, 7.36652525e-07, 9.84413438e-06,
#[Out]#        1.02677275e-04, 8.35896187e-04, 5.31138061e-03, 2.63404079e-02,
#[Out]#        1.01940525e-01, 3.07807651e-01, 7.24955302e-01, 1.33189302e+00,
#[Out]#        1.91005425e+00, 2.14028244e+00, 1.87476215e+00, 1.28304191e+00,
#[Out]#        6.85366883e-01, 2.85575531e-01, 9.28162276e-02, 2.35366198e-02,
#[Out]#        4.65778169e-03, 7.19409735e-04, 8.67263651e-05, 8.16032206e-06,
#[Out]#        5.99301574e-07, 3.43530668e-08, 1.53698236e-09, 5.36667698e-11,
#[Out]#        1.46511902e-12, 3.61757765e-14, 0.00000000e+00, 0.00000000e+00,
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
#[Out]#        9.04392474e-14, 4.37725936e-12, 1.68904322e-10, 5.08767719e-09,
#[Out]#        1.19626013e-07, 2.19567296e-06, 3.14590411e-05, 3.51850993e-04,
#[Out]#        3.07188352e-03, 2.09346826e-02, 1.11344084e-01, 4.61908953e-01,
#[Out]#        1.49268661e+00, 3.75099566e+00, 7.32666065e+00, 1.11594661e+01,
#[Out]#        1.33377158e+01, 1.25648527e+01, 9.31442518e+00, 5.40034875e+00,
#[Out]#        2.43631429e+00, 8.54231656e-01, 2.33108148e-01, 4.95894983e-02,
#[Out]#        8.23129544e-03, 1.06644127e-03, 1.07853261e-04, 8.51460171e-06,
#[Out]#        5.24724236e-07, 2.52425885e-08, 9.47928510e-10, 2.77828939e-11,
#[Out]#        6.33073722e-13, 1.80878197e-14, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        2.71316877e-13, 1.23177856e-11, 4.46732237e-10, 1.26490622e-08,
#[Out]#        2.79605316e-07, 4.82518452e-06, 6.50075414e-05, 6.83745058e-04,
#[Out]#        5.61436058e-03, 3.59876822e-02, 1.80026832e-01, 7.02257606e-01,
#[Out]#        2.13262396e+00, 5.03264483e+00, 9.23182589e+00, 1.32257547e+01,
#[Out]#        1.49030200e+01, 1.32552298e+01, 9.27357304e+00, 5.06728983e+00,
#[Out]#        2.15240233e+00, 7.10444702e-01, 1.82552425e-01, 3.65776485e-02,
#[Out]#        5.71968103e-03, 6.98191868e-04, 6.65354809e-05, 4.95008257e-06,
#[Out]#        2.87509969e-07, 1.30369199e-08, 4.61509317e-10, 1.27518735e-11,
#[Out]#        2.71316444e-13, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        2.17052939e-13, 1.06717690e-11, 4.06087920e-10, 1.20659349e-08,
#[Out]#        2.79916913e-07, 5.07009564e-06, 7.17004250e-05, 7.91671628e-04,
#[Out]#        6.82464388e-03, 4.59292618e-02, 2.41217730e-01, 9.87439226e-01,
#[Out]#        3.14213837e+00, 7.74549928e+00, 1.47824817e+01, 2.19938315e+01,
#[Out]#        2.58318646e+01, 2.41512186e+01, 1.79063291e+01, 1.04003163e+01,
#[Out]#        4.68530445e+00, 1.63364442e+00, 4.42138245e-01, 9.31850910e-02,
#[Out]#        1.53211547e-02, 1.96638793e-03, 1.97038050e-04, 1.54151420e-05,
#[Out]#        9.41593344e-07, 4.49053302e-08, 1.67206481e-09, 4.86016954e-11,
#[Out]#        1.10335068e-12, 1.80877151e-14, 0.00000000e+00, 0.00000000e+00,
#[Out]#        0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
#[Out]#        5.42631265e-14, 3.05682264e-12, 1.35476926e-10, 4.69055821e-09,
#[Out]#        1.26806552e-07, 2.67674112e-06, 4.41182851e-05, 5.67777077e-04,
#[Out]#        5.70530096e-03, 4.47591339e-02, 2.74017724e-01, 1.30661364e+00,
#[Out]#        4.82700837e+00, 1.36814076e+01, 2.94972065e+01, 4.87006446e+01,
#[Out]#        6.34143202e+01, 6.74520501e+01, 5.92848623e+01, 4.21023527e+01,
#[Out]#        2.32926996e+01, 9.81949087e+00, 3.15538419e+00, 7.80856587e-01,
#[Out]#        1.50066428e-01, 2.24838563e-02, 2.62943246e-03, 2.40091940e-04,
#[Out]#        1.71174031e-05, 9.52896076e-07, 4.14209721e-08, 1.48897789e-09,
#[Out]#        2.93671572e-09, 7.88674538e-08, 1.67509034e-06, 2.77812138e-05,
#[Out]#        3.59775940e-04, 3.63812329e-03, 2.87251272e-02, 1.77031478e-01,
#[Out]#        8.50560182e-01, 3.17466589e+00, 9.14474865e+00, 2.02059074e+01,
#[Out]#        3.43789382e+01, 4.59400234e+01, 4.94268927e+01, 4.31970376e+01,
#[Out]#        3.01982462e+01, 1.64619762e+01, 6.89118975e+00, 2.21555061e+00,
#[Out]#        5.50933077e-01, 1.06718613e-01, 1.77700659e-02, 1.67101829e-02,
#[Out]#        1.01509902e-01, 5.40632128e-01, 2.24338737e+00, 7.19495709e+00,
#[Out]#        1.76845147e+01, 3.32590130e+01, 4.86690362e+01, 5.72983288e+01,
#[Out]#        5.63250082e+01, 4.90071031e+01, 4.39471675e+01, 4.76690915e+01,
#[Out]#        5.63346572e+01, 6.06864869e+01, 5.53505801e+01, 4.11952097e+01,
#[Out]#        2.41160954e+01, 1.08297282e+01, 3.71491613e+00, 9.81103683e-01,
#[Out]#        2.01105257e-01, 3.21292172e-02, 4.00658221e-03, 3.90121811e-04,
#[Out]#        2.96624459e-05, 1.76115723e-06, 8.16533718e-08, 2.95620672e-09,
#[Out]#        8.35829523e-11, 1.84493848e-12, 3.61752625e-14, 0.00000000e+00,
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
fitter
#[Out]# <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3cfc084790>
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
ff = 1.2/20 / 5
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3cfc084a00>]
fitter = generate_fitter(modfunc, name="CH3CN")
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=4)
sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=4)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec.quantity))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
ff = 1.2/20 / 5
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d091ff970>]
fitter = generate_fitter(modfunc, name="CH3CN")
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=4)
sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=4)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
sp.plotter()
sp.specfit.plot_fit()
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, np.log10(N_tot.value)])
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
ff = 1.2/20 / 5
ff
#[Out]# 0.012
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec.quantity))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
ff = 1.2/20 / 5
ff
#[Out]# 0.012
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d06bb3fa0>]
fitter = generate_fitter(modfunc, name="CH3CN")
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=4)
sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=4)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, fillingfactor=0.01])
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec.quantity))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
ff = 1.2/20 / 5
ff
#[Out]# 0.012
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d06dd68e0>]
fitter = generate_fitter(modfunc, name="CH3CN")
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
sp.specfit.Registry['CH3CN']
sp.specfit.Registry.multifitters
#[Out]# {'ammonia': <pyspeckit.spectrum.models.ammonia.ammonia_model at 0x2b3d024ebd60>,
#[Out]#  'cold_ammonia': <pyspeckit.spectrum.models.ammonia.cold_ammonia_model at 0x2b3d024ebee0>,
#[Out]#  'ammonia_tau': <pyspeckit.spectrum.models.ammonia.ammonia_model_vtau at 0x2b3d024ebe20>,
#[Out]#  'formaldehyde': <pyspeckit.spectrum.models.formaldehyde.formaldehyde_model at 0x2b3d01a2c340>,
#[Out]#  'gaussian': <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3d024ebf40>,
#[Out]#  'vheightgaussian': <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3d025001f0>,
#[Out]#  'voigt': <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3d02500220>,
#[Out]#  'lorentzian': <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3d025002e0>,
#[Out]#  'ammonia_tau_thin': <pyspeckit.spectrum.models.ammonia.ammonia_model_vtau_thin at 0x2b3d0252c430>,
#[Out]#  'CH3CN': <pyspeckit.spectrum.models.model.SpectralModel at 0x2b3d070cd310>}
sp.specfit.Registry.multifitters['CH3CN'] = fitter
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN']
#[Out]# 5
sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
get_ipython().run_line_magic('debug', '')
fitter.modelfunc
#[Out]# <function __main__.modfunc(xarr, vcen, width, tex, column, fillingfactor)>
get_ipython().run_line_magic('debug', '')
get_ipython().run_line_magic('debug', '')
modfunc
#[Out]# <function __main__.modfunc(xarr, vcen, width, tex, column, fillingfactor)>
get_ipython().run_line_magic('pinfo2', 'modfunc')
get_ipython().run_line_magic('pinfo2', 'generate_fitter')
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
myclass = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule, SpectralModel
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.model SpectralModel
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.model import SpectralModel
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
myclass = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
ff = 1.2/20 / 5
ff
#[Out]# 0.012
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d04426640>]
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"
ff = 1.2/20 / 5
ff
#[Out]# 0.012
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0532d5e0>]
sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01])
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value, 0.01], fixed=[False,False,False,False,True])
sp.plotter()
sp.specfit.plot_fit()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=1):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
ff = 1.2/20 / 5
ff
#[Out]# 0.012
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod * ff)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0376e1c0>]
mod = modfunc(xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d05dba040>]
#sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value], fixed=[False,False,False,False])
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value])
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[v_cen.value, v_disp.value, temp.value, N_tot.value], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =           40  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =            2  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =          167  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =        1e+16  Range:   [0,inf]]
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, N_tot.value], use_lmfit=True)
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =           40  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =            2  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =          167  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =        1e+16  Range:   [0,inf]]
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, N_tot.value], use_lmfit=True, debug=True)
sp.error
#[Out]# masked_array(data=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
#[Out]#                    0., 0., 0., 0., 0., 0., 0., 0.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.error[:]=1
sp.plotter()
sp.specfit.plot_fit()
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, N_tot.value], use_lmfit=True, debug=True)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, N_tot.value], use_lmfit=True, debug=True)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, N_tot.value], use_lmfit=True, debug=True)
sp.error
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.mpfit_status
#[Out]# 'The cosine of the angle between fvec and any column of the jacobian\n  is at most gtol in absolute value.'
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    if N_tot < 100:
        N_tot = 10**N_tot
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, np.log10(N_tot.value)], use_lmfit=True, debug=True)
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    if N_tot < 100:
        N_tot = 10**N_tot
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
ff = 1.2/20 / 5
ff
#[Out]# 0.012
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    if column < 100:
        column = 10**N_tot
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
ff = 1.2/20 / 5
ff
#[Out]# 0.012
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3cfdb13e80>]
#sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.error
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, np.log10(N_tot.value)], use_lmfit=True, debug=True)
get_ipython().run_line_magic('debug', '')
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
ff = 1.2/20 / 5
ff
#[Out]# 0.012
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3cfdbf8e50>]
#sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.error
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit(fittype='CH3CN', guesses=[v_cen.value*1.0, 1.0*v_disp.value, 1.0*temp.value, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit.mpfit_status
#[Out]# 'The cosine of the angle between fvec and any column of the jacobian\n  is at most gtol in absolute value.'
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =           40  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =            2  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =          167  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =           16  Range:   [0,inf]]
sp.error[:]=1
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =         40.1  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =          1.5  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =        155.1  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =           16  Range:   [0,inf]]
sp.error[:]=1
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.residuals
#[Out]# masked_array(data=[-2.19885876e-03, -2.92127689e-03,  6.17697984e-04,
#[Out]#                    -1.04891387e-02, -3.63871907e-02,  5.42679130e-03,
#[Out]#                     5.35453490e-02,  5.00763522e-02,  4.32184920e-02,
#[Out]#                     5.43342921e-02,  6.06385679e-02,  1.86633597e-02,
#[Out]#                    -5.53151470e-02, -6.92727007e-02, -1.81301978e-02,
#[Out]#                     9.93790285e-03, -6.45244242e-03, -1.94484766e-02,
#[Out]#                    -6.73677005e-04,  3.73807317e-02,  5.00296497e-02,
#[Out]#                     2.16284635e-03, -2.58284158e-02, -6.47252514e-03,
#[Out]#                    -2.42833324e-02, -4.75086580e-02, -3.35312650e-02,
#[Out]#                    -2.59206714e-02, -1.21327411e-02, -3.76900423e-03,
#[Out]#                    -8.91143196e-03, -1.10153359e-02, -4.24371627e-02,
#[Out]#                    -7.50772071e-02, -5.39751044e-02, -2.19686038e-02,
#[Out]#                    -3.49981963e-02, -4.78289015e-02, -2.95164084e-02,
#[Out]#                    -9.43871107e-03, -3.38537385e-03, -2.10633296e-02,
#[Out]#                    -4.01293358e-02, -6.17709259e-02, -6.43837305e-02,
#[Out]#                    -4.72865399e-02, -3.97328980e-02, -4.35321825e-02,
#[Out]#                    -3.25133807e-02, -2.44061168e-02, -7.54864836e-02,
#[Out]#                    -8.26471813e-02, -2.85261100e-02, -1.72977914e-02,
#[Out]#                    -2.30501787e-02, -1.36946678e-02,  4.07760191e-03,
#[Out]#                     8.40027402e-03, -2.94669480e-02, -7.48493520e-02,
#[Out]#                    -8.49825609e-02, -6.14444654e-02, -3.48732933e-02,
#[Out]#                    -1.32330924e-02,  4.49348245e-03, -3.50506161e-05,
#[Out]#                    -2.56652265e-02, -3.62311208e-02, -2.01493629e-02,
#[Out]#                     1.11109770e-02,  1.71088626e-02, -2.36234412e-02,
#[Out]#                    -3.93534810e-02, -1.76398052e-02,  1.62198328e-02,
#[Out]#                     3.87501635e-02,  3.35233126e-02,  2.17173027e-02,
#[Out]#                     2.26897186e-02, -1.27061046e-02, -5.94687115e-02,
#[Out]#                    -3.69638954e-02,  1.65152817e-03,  2.53809207e-02,
#[Out]#                     2.20984656e-02, -3.41185002e-02, -8.28873374e-02,
#[Out]#                    -6.01407834e-02, -6.17750461e-03,  4.56074717e-03,
#[Out]#                    -9.53877919e-03, -2.53631020e-02, -3.33464390e-02,
#[Out]#                    -2.02785230e-02, -1.13682623e-02, -2.81508584e-02,
#[Out]#                    -4.80262062e-02, -2.88591844e-02,  5.68329725e-04,
#[Out]#                     1.18813463e-02,  2.75350458e-02,  3.54723666e-02,
#[Out]#                     2.19628459e-02,  1.56750807e-02,  1.10785640e-02,
#[Out]#                    -1.79882233e-02, -5.45882221e-02, -5.96652844e-02,
#[Out]#                    -6.38671954e-02, -8.89373443e-02, -7.21273320e-02,
#[Out]#                    -1.86331792e-02,  3.00853070e-02,  4.56258796e-02,
#[Out]#                     4.76759328e-02,  2.93567797e-02, -1.02592693e-02,
#[Out]#                     2.28277537e-03, -1.02237726e-03, -5.23562846e-02,
#[Out]#                    -5.77521761e-02, -3.68709721e-02, -5.63796595e-02,
#[Out]#                    -5.77448999e-02, -2.29458840e-02, -1.63177121e-02,
#[Out]#                    -2.93555719e-02, -5.04361803e-02, -7.52682878e-02,
#[Out]#                    -6.70144242e-02, -4.51308339e-02, -2.47514046e-02,
#[Out]#                    -1.97764446e-02, -2.62432481e-02, -4.28513767e-02,
#[Out]#                    -3.50528974e-02, -2.79184202e-02, -7.19895750e-02,
#[Out]#                    -7.29488342e-02,  9.30016566e-03,  5.39269429e-02,
#[Out]#                     4.86655968e-02,  3.28674028e-02,  2.48226877e-02,
#[Out]#                     3.75024092e-02,  4.47620348e-02,  2.97478067e-02,
#[Out]#                     2.58450517e-02,  4.07092395e-02,  6.15542252e-02,
#[Out]#                     5.45861628e-02, -6.41937839e-03, -3.36566348e-02,
#[Out]#                     6.64186674e-03,  2.14009994e-02,  1.93306861e-02,
#[Out]#                     2.47626700e-02,  2.04584691e-02,  2.00735220e-02,
#[Out]#                     2.07186286e-02,  1.96046210e-02,  2.63310123e-03,
#[Out]#                    -1.23634495e-02, -2.27966081e-03,  2.21294603e-02,
#[Out]#                     3.15665910e-02,  6.04520969e-02,  8.71186323e-02,
#[Out]#                     4.81207972e-02,  2.86987060e-03, -1.08834338e-02,
#[Out]#                    -3.18026782e-02, -6.09965728e-02, -9.55933358e-02,
#[Out]#                    -1.19476789e-01, -1.08947462e-01, -5.42606062e-02,
#[Out]#                     3.00772658e-02,  6.92333659e-02,  4.45482236e-02,
#[Out]#                     2.64293100e-02,  1.25816635e-02,  1.43552804e-02,
#[Out]#                     7.77722283e-03, -6.50688765e-02, -9.83271054e-02,
#[Out]#                    -3.71520327e-02,  1.64001923e-02,  5.80603426e-04,
#[Out]#                    -3.19288242e-02, -1.19824052e-02,  2.04328778e-02,
#[Out]#                     2.31612461e-02, -4.27018239e-03, -4.07526686e-02,
#[Out]#                    -4.19921144e-02, -2.64471874e-02, -1.06309957e-02,
#[Out]#                    -2.45973051e-03, -1.16185200e-03,  1.09800707e-02,
#[Out]#                     5.91821451e-03, -1.78659852e-02, -1.10376667e-02,
#[Out]#                    -2.47097703e-02, -7.44300901e-02, -6.37065895e-02,
#[Out]#                    -1.86842753e-02, -2.33724919e-02, -3.32966373e-02,
#[Out]#                    -6.03619323e-03,  2.24367883e-02,  1.85465614e-02,
#[Out]#                    -1.34667253e-02, -2.22163647e-02, -1.56370327e-02,
#[Out]#                    -1.24147085e-02, -1.17120418e-02, -2.65478149e-02,
#[Out]#                    -3.20574587e-02, -4.09194769e-03, -1.05930801e-02,
#[Out]#                    -6.05840398e-02, -1.00662190e-01, -1.47102623e-01,
#[Out]#                    -1.81678339e-01, -1.69582786e-01, -1.24273373e-01,
#[Out]#                    -5.39864261e-02,  1.69064494e-02,  2.42851558e-02,
#[Out]#                     5.42846250e-04,  4.96818933e-04,  1.37740957e-03,
#[Out]#                     7.94158818e-03, -1.79799343e-03, -3.22733292e-02,
#[Out]#                    -2.08256037e-02, -9.93382613e-03, -2.35878139e-02,
#[Out]#                    -1.11305782e-02,  2.71435199e-03,  2.79612801e-03,
#[Out]#                    -1.85446580e-02, -6.05544017e-02, -3.21429681e-02,
#[Out]#                     1.40036055e-02, -3.27145070e-02, -6.00152161e-02,
#[Out]#                    -4.52171043e-03,  4.08753925e-02,  5.90340956e-02,
#[Out]#                     5.42433372e-02,  7.12497567e-03, -4.80021851e-02,
#[Out]#                    -4.69873734e-02,  2.58007591e-02,  7.25120086e-02,
#[Out]#                     3.59893764e-02,  5.52693839e-03,  1.08622592e-03,
#[Out]#                    -1.38297713e-03, -3.43643047e-02, -1.21921694e-01,
#[Out]#                    -1.98879795e-01, -2.62305290e-01, -3.36349574e-01,
#[Out]#                    -3.51146994e-01, -2.06271633e-01, -4.62593747e-03,
#[Out]#                     9.03330656e-02,  1.07169381e-01,  6.75664855e-02,
#[Out]#                     3.46292042e-02,  3.19049948e-02,  1.03115223e-02,
#[Out]#                    -4.15994036e-03,  3.13765767e-02,  7.49153596e-02,
#[Out]#                     6.13101364e-02,  2.41080138e-02, -1.64903353e-02,
#[Out]#                    -2.46511041e-02, -2.59439972e-03, -1.31875818e-02,
#[Out]#                    -7.82961846e-03,  6.47090491e-03, -4.10835995e-03,
#[Out]#                     2.31953441e-02,  5.76502419e-02,  2.87194209e-02,
#[Out]#                    -5.50833960e-03, -2.17210343e-02, -7.61297052e-02,
#[Out]#                    -1.80930983e-01, -2.52078327e-01, -2.46741304e-01,
#[Out]#                    -1.79789983e-01, -7.40848476e-02,  2.66898984e-02,
#[Out]#                     9.85795113e-02,  1.24194578e-01,  1.37798325e-01,
#[Out]#                     1.96020972e-01,  2.19266210e-01,  1.84874667e-01,
#[Out]#                     1.63856984e-01,  1.83733526e-01,  2.27764966e-01,
#[Out]#                     2.92536423e-01,  3.33246760e-01,  2.24343645e-01,
#[Out]#                     6.48401883e-03, -1.45709655e-01, -1.27469957e-01,
#[Out]#                     1.28728261e-01,  3.50965215e-01,  2.30161134e-01,
#[Out]#                    -6.82670714e-02, -2.37008985e-01, -1.78938791e-01,
#[Out]#                    -1.09456269e-02,  1.16431726e-01,  1.39140650e-01,
#[Out]#                     6.70012768e-02, -2.38512806e-02, -1.07281358e-01,
#[Out]#                    -1.13477548e-01, -4.20888593e-02,  2.51104694e-02,
#[Out]#                     2.12657477e-03, -7.36871590e-02, -1.10722709e-01,
#[Out]#                    -1.05638975e-01, -7.82027267e-02, -8.11842692e-02,
#[Out]#                    -1.21250954e-01, -1.60518982e-01, -1.64925016e-01,
#[Out]#                    -1.93375811e-01, -1.79895874e-01, -1.40473299e-01,
#[Out]#                    -1.82198772e-01, -2.31527814e-01, -2.42553384e-01,
#[Out]#                    -2.36575042e-01, -2.15408655e-01, -2.21288074e-01,
#[Out]#                    -2.31013052e-01, -1.67138259e-01, -1.20913064e-01,
#[Out]#                    -1.49128586e-01, -1.64463633e-01, -1.48821234e-01,
#[Out]#                    -1.35570624e-01, -1.03112944e-01, -5.62081855e-02,
#[Out]#                    -2.74106900e-02, -1.76595620e-02, -1.32953776e-03,
#[Out]#                     1.91624472e-02,  3.82102631e-03, -1.80389762e-02,
#[Out]#                     1.42930187e-02,  2.60991689e-02,  1.96577587e-02,
#[Out]#                     1.38196876e-02,  1.07092287e-02, -6.66297097e-04,
#[Out]#                    -3.57350550e-02, -3.68685077e-02,  3.50506161e-05,
#[Out]#                     1.28822073e-02,  1.40000467e-02,  2.82818416e-02,
#[Out]#                     7.33185426e-03, -1.00082233e-03,  4.71540733e-02,
#[Out]#                     8.51535640e-02,  6.87381550e-02,  2.78098958e-02,
#[Out]#                     1.87609660e-04, -1.27230770e-02, -5.61549340e-02,
#[Out]#                    -1.04648111e-01, -6.92108926e-02,  5.74925952e-03,
#[Out]#                     1.74143753e-02,  1.36452996e-02,  3.56027056e-02,
#[Out]#                     2.69264718e-02,  1.56107302e-02,  1.78943451e-02,
#[Out]#                     3.44490733e-02,  2.60301866e-02, -1.36545698e-02,
#[Out]#                    -4.20668463e-02, -4.78456547e-02, -6.15690822e-02,
#[Out]#                    -8.08164175e-02],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit.plot_residuals()
sp.specfit.plotresiduals()
sp.xarr
#[Out]# SpectroscopicAxis([147000121483.02917,...,147194951363.59015], unit=Unit("Hz"), refX=None, refX_unit=None, frame=None, redshift=None, xtype=None, velocity convention=None)
sp.xarr.convert_to_unit(u.km/u.s)
freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
sp.xarr.convert_to_unit(u.km/u.s, rest_value=freqs[-1], velocity_convention='radio')
get_ipython().run_line_magic('pinfo', 'sp.xarr.convert_to_unit')
get_ipython().run_line_magic('pinfo', 'sp.xarr.as_unit')
sp.xarr.convert_to_unit_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
sp.xarr
#[Out]# SpectroscopicAxis([355.386323843515,...,-41.47925913444564], unit=Unit("km / s"), refX=<Quantity 147174.5883 MHz>, refX_unit=Unit("MHz"), frame=None, redshift=None, xtype=None, velocity convention=None)
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3cfd818c10>]
#sp.specfit.register_fitter(name='CH3CN', function=fitter, npars=5)
#sp.Registry.add_fitter(name='CH3CN', function=fitter, npars=5)
[v_cen.value, v_disp.value, temp.value, N_tot.value]
#[Out]# [40.0, 2.0, 167.0, 1e+16]
sp.error
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =         40.1  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =          1.5  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =        155.1  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =           16  Range:   [0,inf]]
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =         40.1  Range:[-inf,inf],
#[Out]#  Param #1       WIDTH0 =          1.5  Range:   [0,inf],
#[Out]#  Param #2         TEX0 =        155.1  Range:   [0,inf],
#[Out]#  Param #3      COLUMN0 =           16  Range:   [0,inf]]
sp.error[:]=1
sp.plotter()
sp.specfit.plot_fit()
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], use_lmfit=True, debug=True)
sp.specfit.chi2
#[Out]# 2.8094982929954173
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, np.log10(N_tot.value)], debug=True)
sp.specfit(fittype='CH3CN', guesses=[400.1, 1.5, 15.1, np.log10(N_tot.value)], debug=True)
sp.specfit(fittype='CH3CN', guesses=[400.1, 1.5, 15.1, N_tot.value], debug=True)
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, N_tot.value], debug=True)
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e17], debug=True)
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.data
#[Out]# masked_array(data=[-2.19885876e-03, -2.92127689e-03,  6.17697984e-04,
#[Out]#                    -1.04891387e-02, -3.63871907e-02,  5.42679130e-03,
#[Out]#                     5.35453490e-02,  5.00763522e-02,  4.32184920e-02,
#[Out]#                     5.43342921e-02,  6.06385679e-02,  1.86633597e-02,
#[Out]#                    -5.53151470e-02, -6.92727007e-02, -1.81301978e-02,
#[Out]#                     9.93790285e-03, -6.45244242e-03, -1.94484766e-02,
#[Out]#                    -6.73677005e-04,  3.73807317e-02,  5.00296497e-02,
#[Out]#                     2.16284678e-03, -2.58284085e-02, -6.47242670e-03,
#[Out]#                    -2.42823057e-02, -4.75002991e-02, -3.34781512e-02,
#[Out]#                    -2.56572673e-02, -1.11133358e-02, -6.90927712e-04,
#[Out]#                    -1.66187894e-03,  2.30359425e-03, -2.33366202e-02,
#[Out]#                    -5.36743826e-02, -3.52274829e-02, -9.13818461e-03,
#[Out]#                    -2.81445275e-02, -4.49731462e-02, -2.85882461e-02,
#[Out]#                    -9.20334487e-03, -3.33879603e-03, -2.10561355e-02,
#[Out]#                    -4.01284686e-02, -6.17708443e-02, -6.43837245e-02,
#[Out]#                    -4.72865396e-02, -3.97328980e-02, -4.35321825e-02,
#[Out]#                    -3.25133807e-02, -2.44061168e-02, -7.54864836e-02,
#[Out]#                    -8.26471813e-02, -2.85261100e-02, -1.72977914e-02,
#[Out]#                    -2.30501787e-02, -1.36946678e-02,  4.07760191e-03,
#[Out]#                     8.40027402e-03, -2.94669480e-02, -7.48493520e-02,
#[Out]#                    -8.49825609e-02, -6.14444654e-02, -3.48732933e-02,
#[Out]#                    -1.32330924e-02,  4.49348245e-03, -3.50506161e-05,
#[Out]#                    -2.56652265e-02, -3.62311208e-02, -2.01493629e-02,
#[Out]#                     1.11109770e-02,  1.71088626e-02, -2.36234412e-02,
#[Out]#                    -3.93534810e-02, -1.76398052e-02,  1.62198328e-02,
#[Out]#                     3.87501635e-02,  3.35233126e-02,  2.17173027e-02,
#[Out]#                     2.26897186e-02, -1.27061046e-02, -5.94687115e-02,
#[Out]#                    -3.69638954e-02,  1.65152817e-03,  2.53809207e-02,
#[Out]#                     2.20984656e-02, -3.41185002e-02, -8.28873374e-02,
#[Out]#                    -6.01407834e-02, -6.17750461e-03,  4.56074717e-03,
#[Out]#                    -9.53877919e-03, -2.53631020e-02, -3.33464390e-02,
#[Out]#                    -2.02785230e-02, -1.13682623e-02, -2.81508583e-02,
#[Out]#                    -4.80262050e-02, -2.88591624e-02,  5.68644315e-04,
#[Out]#                     1.18848648e-02,  2.75657646e-02,  3.56817134e-02,
#[Out]#                     2.30762868e-02,  2.02941702e-02,  2.60054301e-02,
#[Out]#                     1.95217334e-02,  1.86783843e-02,  5.19293771e-02,
#[Out]#                     6.95099626e-02,  3.67111822e-02,  2.10169197e-02,
#[Out]#                     3.53703083e-02,  5.44484500e-02,  5.41681962e-02,
#[Out]#                     5.00070143e-02,  2.98526747e-02, -1.01769563e-02,
#[Out]#                     2.29343978e-03, -1.02129873e-03, -5.23561995e-02,
#[Out]#                    -5.77521709e-02, -3.68709719e-02, -5.63796595e-02,
#[Out]#                    -5.77448999e-02, -2.29458840e-02, -1.63177121e-02,
#[Out]#                    -2.93555719e-02, -5.04361803e-02, -7.52682878e-02,
#[Out]#                    -6.70144242e-02, -4.51308339e-02, -2.47514046e-02,
#[Out]#                    -1.97764446e-02, -2.62432481e-02, -4.28513767e-02,
#[Out]#                    -3.50528974e-02, -2.79184202e-02, -7.19895750e-02,
#[Out]#                    -7.29488342e-02,  9.30016566e-03,  5.39269429e-02,
#[Out]#                     4.86655968e-02,  3.28674028e-02,  2.48226877e-02,
#[Out]#                     3.75024092e-02,  4.47620348e-02,  2.97478067e-02,
#[Out]#                     2.58450517e-02,  4.07092395e-02,  6.15542252e-02,
#[Out]#                     5.45861628e-02, -6.41937839e-03, -3.36566348e-02,
#[Out]#                     6.64186674e-03,  2.14009994e-02,  1.93306861e-02,
#[Out]#                     2.47626700e-02,  2.04584691e-02,  2.00735220e-02,
#[Out]#                     2.07186287e-02,  1.96046238e-02,  2.63314948e-03,
#[Out]#                    -1.23627994e-02, -2.27282336e-03,  2.21856039e-02,
#[Out]#                     3.19264678e-02,  6.22523652e-02,  9.41412083e-02,
#[Out]#                     6.94470368e-02,  5.31963189e-02,  8.14348251e-02,
#[Out]#                     1.00454869e-01,  8.80336274e-02,  3.69589620e-02,
#[Out]#                    -2.67410590e-02, -5.82745636e-02, -3.27365829e-02,
#[Out]#                     3.71817128e-02,  7.10588902e-02,  4.49140001e-02,
#[Out]#                     2.64865068e-02,  1.25886454e-02,  1.43559458e-02,
#[Out]#                     7.77727233e-03, -6.50688736e-02, -9.83271053e-02,
#[Out]#                    -3.71520327e-02,  1.64001923e-02,  5.80603426e-04,
#[Out]#                    -3.19288242e-02, -1.19824052e-02,  2.04328778e-02,
#[Out]#                     2.31612461e-02, -4.27018239e-03, -4.07526686e-02,
#[Out]#                    -4.19921144e-02, -2.64471874e-02, -1.06309957e-02,
#[Out]#                    -2.45973051e-03, -1.16185200e-03,  1.09800707e-02,
#[Out]#                     5.91821451e-03, -1.78659852e-02, -1.10376667e-02,
#[Out]#                    -2.47097703e-02, -7.44300901e-02, -6.37065895e-02,
#[Out]#                    -1.86842753e-02, -2.33724919e-02, -3.32966373e-02,
#[Out]#                    -6.03619322e-03,  2.24367884e-02,  1.85465642e-02,
#[Out]#                    -1.34666746e-02, -2.22156477e-02, -1.56291160e-02,
#[Out]#                    -1.23464621e-02, -1.12527492e-02, -2.41356376e-02,
#[Out]#                    -2.21830664e-02,  2.73294360e-02,  6.68619126e-02,
#[Out]#                     8.72407776e-02,  1.19276126e-01,  1.11216023e-01,
#[Out]#                     5.98338463e-02,  9.48050513e-03, -2.02702105e-02,
#[Out]#                    -7.13338156e-03,  3.32428936e-02,  2.87065383e-02,
#[Out]#                     1.47469716e-03,  6.50030481e-04,  1.39707345e-03,
#[Out]#                     7.94355856e-03, -1.79783927e-03, -3.22733198e-02,
#[Out]#                    -2.08256032e-02, -9.93382611e-03, -2.35878139e-02,
#[Out]#                    -1.11305782e-02,  2.71435199e-03,  2.79612801e-03,
#[Out]#                    -1.85446580e-02, -6.05544017e-02, -3.21429681e-02,
#[Out]#                     1.40036055e-02, -3.27145070e-02, -6.00152161e-02,
#[Out]#                    -4.52171043e-03,  4.08753925e-02,  5.90340956e-02,
#[Out]#                     5.42433385e-02,  7.12500244e-03, -4.80017439e-02,
#[Out]#                    -4.69816956e-02,  2.58578121e-02,  7.29595999e-02,
#[Out]#                     3.87295537e-02,  1.85930748e-02,  4.93563096e-02,
#[Out]#                     1.35431098e-01,  2.60607761e-01,  3.65084753e-01,
#[Out]#                     4.35263407e-01,  4.12215212e-01,  2.56499048e-01,
#[Out]#                     6.98765327e-02,  2.66553628e-02,  9.35689712e-02,
#[Out]#                     1.21886908e-01,  1.14977947e-01,  6.90671498e-02,
#[Out]#                     3.48540428e-02,  3.19312892e-02,  1.03139232e-02,
#[Out]#                    -4.15976919e-03,  3.13765863e-02,  7.49153600e-02,
#[Out]#                     6.13101364e-02,  2.41080138e-02, -1.64903345e-02,
#[Out]#                    -2.46510873e-02, -2.59412190e-03, -1.31839840e-02,
#[Out]#                    -7.79323722e-03,  6.75815619e-03, -2.33804517e-03,
#[Out]#                     3.17009459e-02,  8.93969008e-02,  1.20166907e-01,
#[Out]#                     1.96550734e-01,  3.22068347e-01,  3.83270529e-01,
#[Out]#                     3.13337944e-01,  1.79892049e-01,  5.52411586e-02,
#[Out]#                    -1.51702211e-02, -5.17295013e-03,  4.88454046e-02,
#[Out]#                     1.04088842e-01,  1.25261764e-01,  1.37976026e-01,
#[Out]#                     1.96188074e-01,  2.20281309e-01,  1.90280988e-01,
#[Out]#                     1.86290858e-01,  2.55683097e-01,  4.04610113e-01,
#[Out]#                     6.25126554e-01,  8.19937122e-01,  7.97326933e-01,
#[Out]#                     5.69734101e-01,  3.44361376e-01,  3.12001717e-01,
#[Out]#                     6.05419176e-01,  9.14311787e-01,  8.37026004e-01,
#[Out]#                     4.85238729e-01,  1.74943112e-01,  6.22221633e-02,
#[Out]#                     9.73516552e-02,  1.53580887e-01,  1.48951686e-01,
#[Out]#                     6.90123293e-02, -2.35299885e-02, -1.07241292e-01,
#[Out]#                    -1.13473646e-01, -4.20885626e-02,  2.51104870e-02,
#[Out]#                     2.12657559e-03, -7.36871590e-02, -1.10722709e-01,
#[Out]#                    -1.05638975e-01, -7.82027267e-02, -8.11842692e-02,
#[Out]#                    -1.21250954e-01, -1.60518982e-01, -1.64925016e-01,
#[Out]#                    -1.93375811e-01, -1.79895874e-01, -1.40473299e-01,
#[Out]#                    -1.82198772e-01, -2.31527814e-01, -2.42553384e-01,
#[Out]#                    -2.36575042e-01, -2.15408655e-01, -2.21288074e-01,
#[Out]#                    -2.31013052e-01, -1.67138259e-01, -1.20913064e-01,
#[Out]#                    -1.49128586e-01, -1.64463633e-01, -1.48821234e-01,
#[Out]#                    -1.35570624e-01, -1.03112944e-01, -5.62081855e-02,
#[Out]#                    -2.74106900e-02, -1.76595620e-02, -1.32953776e-03,
#[Out]#                     1.91624472e-02,  3.82102631e-03, -1.80389762e-02,
#[Out]#                     1.42930187e-02,  2.60991689e-02,  1.96577587e-02,
#[Out]#                     1.38196876e-02,  1.07092287e-02, -6.66297097e-04,
#[Out]#                    -3.57350550e-02, -3.68685077e-02,  3.50506161e-05,
#[Out]#                     1.28822073e-02,  1.40000467e-02,  2.82818416e-02,
#[Out]#                     7.33185426e-03, -1.00082233e-03,  4.71540733e-02,
#[Out]#                     8.51535640e-02,  6.87381550e-02,  2.78098958e-02,
#[Out]#                     1.87609660e-04, -1.27230770e-02, -5.61549340e-02,
#[Out]#                    -1.04648111e-01, -6.92108926e-02,  5.74925952e-03,
#[Out]#                     1.74143753e-02,  1.36452996e-02,  3.56027056e-02,
#[Out]#                     2.69264718e-02,  1.56107302e-02,  1.78943451e-02,
#[Out]#                     3.44490733e-02,  2.60301866e-02, -1.36545698e-02,
#[Out]#                    -4.20668463e-02, -4.78456547e-02, -6.15690822e-02,
#[Out]#                    -8.08164175e-02],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.xmin,sp.xmax
sp.specfit.xmin,sp.specfit.xmax
#[Out]# (0, 400)
sp.data[0:400]
#[Out]# masked_array(data=[-2.19885876e-03, -2.92127689e-03,  6.17697984e-04,
#[Out]#                    -1.04891387e-02, -3.63871907e-02,  5.42679130e-03,
#[Out]#                     5.35453490e-02,  5.00763522e-02,  4.32184920e-02,
#[Out]#                     5.43342921e-02,  6.06385679e-02,  1.86633597e-02,
#[Out]#                    -5.53151470e-02, -6.92727007e-02, -1.81301978e-02,
#[Out]#                     9.93790285e-03, -6.45244242e-03, -1.94484766e-02,
#[Out]#                    -6.73677005e-04,  3.73807317e-02,  5.00296497e-02,
#[Out]#                     2.16284678e-03, -2.58284085e-02, -6.47242670e-03,
#[Out]#                    -2.42823057e-02, -4.75002991e-02, -3.34781512e-02,
#[Out]#                    -2.56572673e-02, -1.11133358e-02, -6.90927712e-04,
#[Out]#                    -1.66187894e-03,  2.30359425e-03, -2.33366202e-02,
#[Out]#                    -5.36743826e-02, -3.52274829e-02, -9.13818461e-03,
#[Out]#                    -2.81445275e-02, -4.49731462e-02, -2.85882461e-02,
#[Out]#                    -9.20334487e-03, -3.33879603e-03, -2.10561355e-02,
#[Out]#                    -4.01284686e-02, -6.17708443e-02, -6.43837245e-02,
#[Out]#                    -4.72865396e-02, -3.97328980e-02, -4.35321825e-02,
#[Out]#                    -3.25133807e-02, -2.44061168e-02, -7.54864836e-02,
#[Out]#                    -8.26471813e-02, -2.85261100e-02, -1.72977914e-02,
#[Out]#                    -2.30501787e-02, -1.36946678e-02,  4.07760191e-03,
#[Out]#                     8.40027402e-03, -2.94669480e-02, -7.48493520e-02,
#[Out]#                    -8.49825609e-02, -6.14444654e-02, -3.48732933e-02,
#[Out]#                    -1.32330924e-02,  4.49348245e-03, -3.50506161e-05,
#[Out]#                    -2.56652265e-02, -3.62311208e-02, -2.01493629e-02,
#[Out]#                     1.11109770e-02,  1.71088626e-02, -2.36234412e-02,
#[Out]#                    -3.93534810e-02, -1.76398052e-02,  1.62198328e-02,
#[Out]#                     3.87501635e-02,  3.35233126e-02,  2.17173027e-02,
#[Out]#                     2.26897186e-02, -1.27061046e-02, -5.94687115e-02,
#[Out]#                    -3.69638954e-02,  1.65152817e-03,  2.53809207e-02,
#[Out]#                     2.20984656e-02, -3.41185002e-02, -8.28873374e-02,
#[Out]#                    -6.01407834e-02, -6.17750461e-03,  4.56074717e-03,
#[Out]#                    -9.53877919e-03, -2.53631020e-02, -3.33464390e-02,
#[Out]#                    -2.02785230e-02, -1.13682623e-02, -2.81508583e-02,
#[Out]#                    -4.80262050e-02, -2.88591624e-02,  5.68644315e-04,
#[Out]#                     1.18848648e-02,  2.75657646e-02,  3.56817134e-02,
#[Out]#                     2.30762868e-02,  2.02941702e-02,  2.60054301e-02,
#[Out]#                     1.95217334e-02,  1.86783843e-02,  5.19293771e-02,
#[Out]#                     6.95099626e-02,  3.67111822e-02,  2.10169197e-02,
#[Out]#                     3.53703083e-02,  5.44484500e-02,  5.41681962e-02,
#[Out]#                     5.00070143e-02,  2.98526747e-02, -1.01769563e-02,
#[Out]#                     2.29343978e-03, -1.02129873e-03, -5.23561995e-02,
#[Out]#                    -5.77521709e-02, -3.68709719e-02, -5.63796595e-02,
#[Out]#                    -5.77448999e-02, -2.29458840e-02, -1.63177121e-02,
#[Out]#                    -2.93555719e-02, -5.04361803e-02, -7.52682878e-02,
#[Out]#                    -6.70144242e-02, -4.51308339e-02, -2.47514046e-02,
#[Out]#                    -1.97764446e-02, -2.62432481e-02, -4.28513767e-02,
#[Out]#                    -3.50528974e-02, -2.79184202e-02, -7.19895750e-02,
#[Out]#                    -7.29488342e-02,  9.30016566e-03,  5.39269429e-02,
#[Out]#                     4.86655968e-02,  3.28674028e-02,  2.48226877e-02,
#[Out]#                     3.75024092e-02,  4.47620348e-02,  2.97478067e-02,
#[Out]#                     2.58450517e-02,  4.07092395e-02,  6.15542252e-02,
#[Out]#                     5.45861628e-02, -6.41937839e-03, -3.36566348e-02,
#[Out]#                     6.64186674e-03,  2.14009994e-02,  1.93306861e-02,
#[Out]#                     2.47626700e-02,  2.04584691e-02,  2.00735220e-02,
#[Out]#                     2.07186287e-02,  1.96046238e-02,  2.63314948e-03,
#[Out]#                    -1.23627994e-02, -2.27282336e-03,  2.21856039e-02,
#[Out]#                     3.19264678e-02,  6.22523652e-02,  9.41412083e-02,
#[Out]#                     6.94470368e-02,  5.31963189e-02,  8.14348251e-02,
#[Out]#                     1.00454869e-01,  8.80336274e-02,  3.69589620e-02,
#[Out]#                    -2.67410590e-02, -5.82745636e-02, -3.27365829e-02,
#[Out]#                     3.71817128e-02,  7.10588902e-02,  4.49140001e-02,
#[Out]#                     2.64865068e-02,  1.25886454e-02,  1.43559458e-02,
#[Out]#                     7.77727233e-03, -6.50688736e-02, -9.83271053e-02,
#[Out]#                    -3.71520327e-02,  1.64001923e-02,  5.80603426e-04,
#[Out]#                    -3.19288242e-02, -1.19824052e-02,  2.04328778e-02,
#[Out]#                     2.31612461e-02, -4.27018239e-03, -4.07526686e-02,
#[Out]#                    -4.19921144e-02, -2.64471874e-02, -1.06309957e-02,
#[Out]#                    -2.45973051e-03, -1.16185200e-03,  1.09800707e-02,
#[Out]#                     5.91821451e-03, -1.78659852e-02, -1.10376667e-02,
#[Out]#                    -2.47097703e-02, -7.44300901e-02, -6.37065895e-02,
#[Out]#                    -1.86842753e-02, -2.33724919e-02, -3.32966373e-02,
#[Out]#                    -6.03619322e-03,  2.24367884e-02,  1.85465642e-02,
#[Out]#                    -1.34666746e-02, -2.22156477e-02, -1.56291160e-02,
#[Out]#                    -1.23464621e-02, -1.12527492e-02, -2.41356376e-02,
#[Out]#                    -2.21830664e-02,  2.73294360e-02,  6.68619126e-02,
#[Out]#                     8.72407776e-02,  1.19276126e-01,  1.11216023e-01,
#[Out]#                     5.98338463e-02,  9.48050513e-03, -2.02702105e-02,
#[Out]#                    -7.13338156e-03,  3.32428936e-02,  2.87065383e-02,
#[Out]#                     1.47469716e-03,  6.50030481e-04,  1.39707345e-03,
#[Out]#                     7.94355856e-03, -1.79783927e-03, -3.22733198e-02,
#[Out]#                    -2.08256032e-02, -9.93382611e-03, -2.35878139e-02,
#[Out]#                    -1.11305782e-02,  2.71435199e-03,  2.79612801e-03,
#[Out]#                    -1.85446580e-02, -6.05544017e-02, -3.21429681e-02,
#[Out]#                     1.40036055e-02, -3.27145070e-02, -6.00152161e-02,
#[Out]#                    -4.52171043e-03,  4.08753925e-02,  5.90340956e-02,
#[Out]#                     5.42433385e-02,  7.12500244e-03, -4.80017439e-02,
#[Out]#                    -4.69816956e-02,  2.58578121e-02,  7.29595999e-02,
#[Out]#                     3.87295537e-02,  1.85930748e-02,  4.93563096e-02,
#[Out]#                     1.35431098e-01,  2.60607761e-01,  3.65084753e-01,
#[Out]#                     4.35263407e-01,  4.12215212e-01,  2.56499048e-01,
#[Out]#                     6.98765327e-02,  2.66553628e-02,  9.35689712e-02,
#[Out]#                     1.21886908e-01,  1.14977947e-01,  6.90671498e-02,
#[Out]#                     3.48540428e-02,  3.19312892e-02,  1.03139232e-02,
#[Out]#                    -4.15976919e-03,  3.13765863e-02,  7.49153600e-02,
#[Out]#                     6.13101364e-02,  2.41080138e-02, -1.64903345e-02,
#[Out]#                    -2.46510873e-02, -2.59412190e-03, -1.31839840e-02,
#[Out]#                    -7.79323722e-03,  6.75815619e-03, -2.33804517e-03,
#[Out]#                     3.17009459e-02,  8.93969008e-02,  1.20166907e-01,
#[Out]#                     1.96550734e-01,  3.22068347e-01,  3.83270529e-01,
#[Out]#                     3.13337944e-01,  1.79892049e-01,  5.52411586e-02,
#[Out]#                    -1.51702211e-02, -5.17295013e-03,  4.88454046e-02,
#[Out]#                     1.04088842e-01,  1.25261764e-01,  1.37976026e-01,
#[Out]#                     1.96188074e-01,  2.20281309e-01,  1.90280988e-01,
#[Out]#                     1.86290858e-01,  2.55683097e-01,  4.04610113e-01,
#[Out]#                     6.25126554e-01,  8.19937122e-01,  7.97326933e-01,
#[Out]#                     5.69734101e-01,  3.44361376e-01,  3.12001717e-01,
#[Out]#                     6.05419176e-01,  9.14311787e-01,  8.37026004e-01,
#[Out]#                     4.85238729e-01,  1.74943112e-01,  6.22221633e-02,
#[Out]#                     9.73516552e-02,  1.53580887e-01,  1.48951686e-01,
#[Out]#                     6.90123293e-02, -2.35299885e-02, -1.07241292e-01,
#[Out]#                    -1.13473646e-01, -4.20885626e-02,  2.51104870e-02,
#[Out]#                     2.12657559e-03, -7.36871590e-02, -1.10722709e-01,
#[Out]#                    -1.05638975e-01, -7.82027267e-02, -8.11842692e-02,
#[Out]#                    -1.21250954e-01, -1.60518982e-01, -1.64925016e-01,
#[Out]#                    -1.93375811e-01, -1.79895874e-01, -1.40473299e-01,
#[Out]#                    -1.82198772e-01, -2.31527814e-01, -2.42553384e-01,
#[Out]#                    -2.36575042e-01, -2.15408655e-01, -2.21288074e-01,
#[Out]#                    -2.31013052e-01, -1.67138259e-01, -1.20913064e-01,
#[Out]#                    -1.49128586e-01, -1.64463633e-01, -1.48821234e-01,
#[Out]#                    -1.35570624e-01, -1.03112944e-01, -5.62081855e-02,
#[Out]#                    -2.74106900e-02, -1.76595620e-02, -1.32953776e-03,
#[Out]#                     1.91624472e-02,  3.82102631e-03, -1.80389762e-02,
#[Out]#                     1.42930187e-02,  2.60991689e-02,  1.96577587e-02,
#[Out]#                     1.38196876e-02,  1.07092287e-02, -6.66297097e-04,
#[Out]#                    -3.57350550e-02, -3.68685077e-02,  3.50506161e-05,
#[Out]#                     1.28822073e-02,  1.40000467e-02,  2.82818416e-02,
#[Out]#                     7.33185426e-03, -1.00082233e-03,  4.71540733e-02,
#[Out]#                     8.51535640e-02,  6.87381550e-02,  2.78098958e-02,
#[Out]#                     1.87609660e-04, -1.27230770e-02, -5.61549340e-02,
#[Out]#                    -1.04648111e-01, -6.92108926e-02,  5.74925952e-03,
#[Out]#                     1.74143753e-02,  1.36452996e-02,  3.56027056e-02,
#[Out]#                     2.69264718e-02,  1.56107302e-02,  1.78943451e-02,
#[Out]#                     3.44490733e-02,  2.60301866e-02, -1.36545698e-02,
#[Out]#                    -4.20668463e-02, -4.78456547e-02, -6.15690822e-02,
#[Out]#                    -8.08164175e-02],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit.spectofit[0:400]
#[Out]# masked_array(data=[-2.19885876e-03, -2.92127689e-03,  6.17697984e-04,
#[Out]#                    -1.04891387e-02, -3.63871907e-02,  5.42679130e-03,
#[Out]#                     5.35453490e-02,  5.00763522e-02,  4.32184920e-02,
#[Out]#                     5.43342921e-02,  6.06385679e-02,  1.86633597e-02,
#[Out]#                    -5.53151470e-02, -6.92727007e-02, -1.81301978e-02,
#[Out]#                     9.93790285e-03, -6.45244242e-03, -1.94484766e-02,
#[Out]#                    -6.73677005e-04,  3.73807317e-02,  5.00296497e-02,
#[Out]#                     2.16284678e-03, -2.58284085e-02, -6.47242670e-03,
#[Out]#                    -2.42823057e-02, -4.75002991e-02, -3.34781512e-02,
#[Out]#                    -2.56572673e-02, -1.11133358e-02, -6.90927712e-04,
#[Out]#                    -1.66187894e-03,  2.30359425e-03, -2.33366202e-02,
#[Out]#                    -5.36743826e-02, -3.52274829e-02, -9.13818461e-03,
#[Out]#                    -2.81445275e-02, -4.49731462e-02, -2.85882461e-02,
#[Out]#                    -9.20334487e-03, -3.33879603e-03, -2.10561355e-02,
#[Out]#                    -4.01284686e-02, -6.17708443e-02, -6.43837245e-02,
#[Out]#                    -4.72865396e-02, -3.97328980e-02, -4.35321825e-02,
#[Out]#                    -3.25133807e-02, -2.44061168e-02, -7.54864836e-02,
#[Out]#                    -8.26471813e-02, -2.85261100e-02, -1.72977914e-02,
#[Out]#                    -2.30501787e-02, -1.36946678e-02,  4.07760191e-03,
#[Out]#                     8.40027402e-03, -2.94669480e-02, -7.48493520e-02,
#[Out]#                    -8.49825609e-02, -6.14444654e-02, -3.48732933e-02,
#[Out]#                    -1.32330924e-02,  4.49348245e-03, -3.50506161e-05,
#[Out]#                    -2.56652265e-02, -3.62311208e-02, -2.01493629e-02,
#[Out]#                     1.11109770e-02,  1.71088626e-02, -2.36234412e-02,
#[Out]#                    -3.93534810e-02, -1.76398052e-02,  1.62198328e-02,
#[Out]#                     3.87501635e-02,  3.35233126e-02,  2.17173027e-02,
#[Out]#                     2.26897186e-02, -1.27061046e-02, -5.94687115e-02,
#[Out]#                    -3.69638954e-02,  1.65152817e-03,  2.53809207e-02,
#[Out]#                     2.20984656e-02, -3.41185002e-02, -8.28873374e-02,
#[Out]#                    -6.01407834e-02, -6.17750461e-03,  4.56074717e-03,
#[Out]#                    -9.53877919e-03, -2.53631020e-02, -3.33464390e-02,
#[Out]#                    -2.02785230e-02, -1.13682623e-02, -2.81508583e-02,
#[Out]#                    -4.80262050e-02, -2.88591624e-02,  5.68644315e-04,
#[Out]#                     1.18848648e-02,  2.75657646e-02,  3.56817134e-02,
#[Out]#                     2.30762868e-02,  2.02941702e-02,  2.60054301e-02,
#[Out]#                     1.95217334e-02,  1.86783843e-02,  5.19293771e-02,
#[Out]#                     6.95099626e-02,  3.67111822e-02,  2.10169197e-02,
#[Out]#                     3.53703083e-02,  5.44484500e-02,  5.41681962e-02,
#[Out]#                     5.00070143e-02,  2.98526747e-02, -1.01769563e-02,
#[Out]#                     2.29343978e-03, -1.02129873e-03, -5.23561995e-02,
#[Out]#                    -5.77521709e-02, -3.68709719e-02, -5.63796595e-02,
#[Out]#                    -5.77448999e-02, -2.29458840e-02, -1.63177121e-02,
#[Out]#                    -2.93555719e-02, -5.04361803e-02, -7.52682878e-02,
#[Out]#                    -6.70144242e-02, -4.51308339e-02, -2.47514046e-02,
#[Out]#                    -1.97764446e-02, -2.62432481e-02, -4.28513767e-02,
#[Out]#                    -3.50528974e-02, -2.79184202e-02, -7.19895750e-02,
#[Out]#                    -7.29488342e-02,  9.30016566e-03,  5.39269429e-02,
#[Out]#                     4.86655968e-02,  3.28674028e-02,  2.48226877e-02,
#[Out]#                     3.75024092e-02,  4.47620348e-02,  2.97478067e-02,
#[Out]#                     2.58450517e-02,  4.07092395e-02,  6.15542252e-02,
#[Out]#                     5.45861628e-02, -6.41937839e-03, -3.36566348e-02,
#[Out]#                     6.64186674e-03,  2.14009994e-02,  1.93306861e-02,
#[Out]#                     2.47626700e-02,  2.04584691e-02,  2.00735220e-02,
#[Out]#                     2.07186287e-02,  1.96046238e-02,  2.63314948e-03,
#[Out]#                    -1.23627994e-02, -2.27282336e-03,  2.21856039e-02,
#[Out]#                     3.19264678e-02,  6.22523652e-02,  9.41412083e-02,
#[Out]#                     6.94470368e-02,  5.31963189e-02,  8.14348251e-02,
#[Out]#                     1.00454869e-01,  8.80336274e-02,  3.69589620e-02,
#[Out]#                    -2.67410590e-02, -5.82745636e-02, -3.27365829e-02,
#[Out]#                     3.71817128e-02,  7.10588902e-02,  4.49140001e-02,
#[Out]#                     2.64865068e-02,  1.25886454e-02,  1.43559458e-02,
#[Out]#                     7.77727233e-03, -6.50688736e-02, -9.83271053e-02,
#[Out]#                    -3.71520327e-02,  1.64001923e-02,  5.80603426e-04,
#[Out]#                    -3.19288242e-02, -1.19824052e-02,  2.04328778e-02,
#[Out]#                     2.31612461e-02, -4.27018239e-03, -4.07526686e-02,
#[Out]#                    -4.19921144e-02, -2.64471874e-02, -1.06309957e-02,
#[Out]#                    -2.45973051e-03, -1.16185200e-03,  1.09800707e-02,
#[Out]#                     5.91821451e-03, -1.78659852e-02, -1.10376667e-02,
#[Out]#                    -2.47097703e-02, -7.44300901e-02, -6.37065895e-02,
#[Out]#                    -1.86842753e-02, -2.33724919e-02, -3.32966373e-02,
#[Out]#                    -6.03619322e-03,  2.24367884e-02,  1.85465642e-02,
#[Out]#                    -1.34666746e-02, -2.22156477e-02, -1.56291160e-02,
#[Out]#                    -1.23464621e-02, -1.12527492e-02, -2.41356376e-02,
#[Out]#                    -2.21830664e-02,  2.73294360e-02,  6.68619126e-02,
#[Out]#                     8.72407776e-02,  1.19276126e-01,  1.11216023e-01,
#[Out]#                     5.98338463e-02,  9.48050513e-03, -2.02702105e-02,
#[Out]#                    -7.13338156e-03,  3.32428936e-02,  2.87065383e-02,
#[Out]#                     1.47469716e-03,  6.50030481e-04,  1.39707345e-03,
#[Out]#                     7.94355856e-03, -1.79783927e-03, -3.22733198e-02,
#[Out]#                    -2.08256032e-02, -9.93382611e-03, -2.35878139e-02,
#[Out]#                    -1.11305782e-02,  2.71435199e-03,  2.79612801e-03,
#[Out]#                    -1.85446580e-02, -6.05544017e-02, -3.21429681e-02,
#[Out]#                     1.40036055e-02, -3.27145070e-02, -6.00152161e-02,
#[Out]#                    -4.52171043e-03,  4.08753925e-02,  5.90340956e-02,
#[Out]#                     5.42433385e-02,  7.12500244e-03, -4.80017439e-02,
#[Out]#                    -4.69816956e-02,  2.58578121e-02,  7.29595999e-02,
#[Out]#                     3.87295537e-02,  1.85930748e-02,  4.93563096e-02,
#[Out]#                     1.35431098e-01,  2.60607761e-01,  3.65084753e-01,
#[Out]#                     4.35263407e-01,  4.12215212e-01,  2.56499048e-01,
#[Out]#                     6.98765327e-02,  2.66553628e-02,  9.35689712e-02,
#[Out]#                     1.21886908e-01,  1.14977947e-01,  6.90671498e-02,
#[Out]#                     3.48540428e-02,  3.19312892e-02,  1.03139232e-02,
#[Out]#                    -4.15976919e-03,  3.13765863e-02,  7.49153600e-02,
#[Out]#                     6.13101364e-02,  2.41080138e-02, -1.64903345e-02,
#[Out]#                    -2.46510873e-02, -2.59412190e-03, -1.31839840e-02,
#[Out]#                    -7.79323722e-03,  6.75815619e-03, -2.33804517e-03,
#[Out]#                     3.17009459e-02,  8.93969008e-02,  1.20166907e-01,
#[Out]#                     1.96550734e-01,  3.22068347e-01,  3.83270529e-01,
#[Out]#                     3.13337944e-01,  1.79892049e-01,  5.52411586e-02,
#[Out]#                    -1.51702211e-02, -5.17295013e-03,  4.88454046e-02,
#[Out]#                     1.04088842e-01,  1.25261764e-01,  1.37976026e-01,
#[Out]#                     1.96188074e-01,  2.20281309e-01,  1.90280988e-01,
#[Out]#                     1.86290858e-01,  2.55683097e-01,  4.04610113e-01,
#[Out]#                     6.25126554e-01,  8.19937122e-01,  7.97326933e-01,
#[Out]#                     5.69734101e-01,  3.44361376e-01,  3.12001717e-01,
#[Out]#                     6.05419176e-01,  9.14311787e-01,  8.37026004e-01,
#[Out]#                     4.85238729e-01,  1.74943112e-01,  6.22221633e-02,
#[Out]#                     9.73516552e-02,  1.53580887e-01,  1.48951686e-01,
#[Out]#                     6.90123293e-02, -2.35299885e-02, -1.07241292e-01,
#[Out]#                    -1.13473646e-01, -4.20885626e-02,  2.51104870e-02,
#[Out]#                     2.12657559e-03, -7.36871590e-02, -1.10722709e-01,
#[Out]#                    -1.05638975e-01, -7.82027267e-02, -8.11842692e-02,
#[Out]#                    -1.21250954e-01, -1.60518982e-01, -1.64925016e-01,
#[Out]#                    -1.93375811e-01, -1.79895874e-01, -1.40473299e-01,
#[Out]#                    -1.82198772e-01, -2.31527814e-01, -2.42553384e-01,
#[Out]#                    -2.36575042e-01, -2.15408655e-01, -2.21288074e-01,
#[Out]#                    -2.31013052e-01, -1.67138259e-01, -1.20913064e-01,
#[Out]#                    -1.49128586e-01, -1.64463633e-01, -1.48821234e-01,
#[Out]#                    -1.35570624e-01, -1.03112944e-01, -5.62081855e-02,
#[Out]#                    -2.74106900e-02, -1.76595620e-02, -1.32953776e-03,
#[Out]#                     1.91624472e-02,  3.82102631e-03, -1.80389762e-02,
#[Out]#                     1.42930187e-02,  2.60991689e-02,  1.96577587e-02,
#[Out]#                     1.38196876e-02,  1.07092287e-02, -6.66297097e-04,
#[Out]#                    -3.57350550e-02, -3.68685077e-02,  3.50506161e-05,
#[Out]#                     1.28822073e-02,  1.40000467e-02,  2.82818416e-02,
#[Out]#                     7.33185426e-03, -1.00082233e-03,  4.71540733e-02,
#[Out]#                     8.51535640e-02,  6.87381550e-02,  2.78098958e-02,
#[Out]#                     1.87609660e-04, -1.27230770e-02, -5.61549340e-02,
#[Out]#                    -1.04648111e-01, -6.92108926e-02,  5.74925952e-03,
#[Out]#                     1.74143753e-02,  1.36452996e-02,  3.56027056e-02,
#[Out]#                     2.69264718e-02,  1.56107302e-02,  1.78943451e-02,
#[Out]#                     3.44490733e-02,  2.60301866e-02, -1.36545698e-02,
#[Out]#                    -4.20668463e-02, -4.78456547e-02, -6.15690822e-02,
#[Out]#                    -8.08164175e-02],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
sp.specfit.errspec
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
from astropy import log
log.setLevel('DEBUG')
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 100, 15))
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 200, 14))
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0bb59b20>]
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 100, 15))
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 200, 1e14))
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0bc4d8b0>]
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.01):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(sp.xarr, vcen, width, tex, column,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 100, 15))
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 200, 1e14))
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0bcff850>]
freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
from astropy import log
log.setLevel('INFO')
sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.specfit.chi2
#[Out]# 1.9024683591211478
sp.specfit.parinfo
#[Out]# [Param #0       SHIFT0 =      40.9587 +/-         1.24298 ,
#[Out]#  Param #1       WIDTH0 =      1.51585 +/-         1.37145   Range:   [0,inf),
#[Out]#  Param #2         TEX0 =      79.9671 +/-         140.495   Range:   [0,inf),
#[Out]#  Param #3      COLUMN0 =  5.07044e+15 +/-               0   Range:   [0,inf)]
sp.error[:]=1
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
sp.specfit.xmin,sp.specfit.xmax
#[Out]# (0, 400)
sp.specfit.errspec
#[Out]# masked_array(data=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
#[Out]#                    1., 1., 1., 1., 1., 1., 1., 1.],
#[Out]#              mask=False,
#[Out]#        fill_value=1e+20)
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15, 0.01])
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e17, 0.01], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e16, 0.01], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15, 0.01], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15, 0.01], use_lmfit=False)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.01], use_lmfit=False)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15, 0.01], use_lmfit=False)
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15, 0.01], use_lmfit=False)
sp.plotter()
sp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.01], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.012):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(sp.xarr, vcen, width, tex, column,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 100, 15))
pl.plot(sp.xarr, modfunc(sp.xarr, 40.1, 1, 200, 1e14))
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d125e5fa0>]
freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
sp.xarr
#[Out]# SpectroscopicAxis([355.386323843515,...,-41.47925913444564], unit=Unit("km / s"), refX=<Quantity 147174.5883 MHz>, refX_unit=Unit("MHz"), frame=None, redshift=None, xtype=None, velocity convention=None)
freqs
#[Out]# <Quantity [147035.8351, 147072.6021, 147103.738 , 147129.2302, 147149.0683,
#[Out]#            147163.2441, 147171.7519, 147174.5883] MHz>
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
sp.xarr
#[Out]# SpectroscopicAxis([355.386323843515,...,-41.47925913444564], unit=Unit("km / s"), refX=<Quantity 147174.5883 MHz>, refX_unit=Unit("MHz"), frame=None, redshift=None, xtype=None, velocity convention=None)
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d12787fa0>]
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d127f3b80>]
fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.012], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.02):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(sp.xarr, vcen, width, tex, column,
                                  freqs, aij, deg, EU, partfunc)*fillingfactor
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d12b9db50>]
fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.012], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
dcube = SpectralCube.read('/orange/adamginsburg/brick_alma_linesurvey/danwalker/Brick_CH3CN.image.fits')
reg = regions.Regions.read('../centralcoreellipse.reg')
dscube = dcube.subcube_from_regions(reg)
dscube
#[Out]# SpectralCube with shape=(700, 117, 114) and unit=Jy / beam:
#[Out]#  n_x:    114  type_x: RA---SIN  unit_x: deg    range:   266.543732 deg:  266.544806 deg
#[Out]#  n_y:    117  type_y: DEC--SIN  unit_y: deg    range:   -28.705359 deg:  -28.704392 deg
#[Out]#  n_s:    700  type_s: VRAD      unit_s: m / s  range:   -63094.684 m / s:  400457.822 m / s
m0 = dscube.spectral_slab(220.71*u.GHz, 220.72*u.GHz).moment0() # slab of comps 0 and 1
m0 = dscube.with_spectral_unit(u.GHz).spectral_slab(220.71*u.GHz, 220.72*u.GHz).moment0() # slab of comps 0 and 1
m0.quicklook()
davgspec = (scube.spectral_slab(220.71*u.GHz, 220.72*u.GHz) * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
davgspec = (dscube.spectral_slab(220.71*u.GHz, 220.72*u.GHz) * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
davgspec = (dscube * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
davgspec
#[Out]# <OneDSpectrum [-1.52124130e-05, 1.37156545e-04, 1.69312334e-04,
#[Out]#                -5.06345132e-05,-1.77297974e-04,-1.02236598e-04,
#[Out]#                -4.05921419e-05,-2.57677370e-04,-8.93292163e-05,
#[Out]#                -8.02707800e-05,-2.50306657e-05,-1.60126307e-04,
#[Out]#                -6.24471359e-05, 1.07153879e-04, 1.69636405e-05,
#[Out]#                -7.10016830e-05,-2.52348655e-05,-2.32760256e-04,
#[Out]#                -2.18755653e-04,-3.29219067e-04,-3.86474421e-04,
#[Out]#                -3.35996992e-05, 1.37284442e-04,-3.02511253e-05,
#[Out]#                -1.58321476e-04,-2.52710684e-04,-1.20148768e-04,
#[Out]#                -6.54856631e-05,-4.60070078e-05,-7.25271675e-05,
#[Out]#                -2.77294166e-04,-1.15548923e-04, 8.01606002e-05,
#[Out]#                -1.60352531e-04,-2.58178887e-04,-3.22304259e-05,
#[Out]#                -4.09425447e-05,-6.32898518e-05, 1.02313017e-04,
#[Out]#                 2.58506552e-05, 4.30196924e-05, 1.63583201e-04,
#[Out]#                 2.06414959e-04, 1.89673447e-04,-4.97743395e-06,
#[Out]#                 1.64669418e-05, 1.10300731e-04, 2.01906369e-04,
#[Out]#                 2.31045065e-04, 4.28017927e-04, 3.25111760e-04,
#[Out]#                 2.10234884e-05,-1.56564220e-05, 3.07351729e-04,
#[Out]#                 2.08330763e-04, 2.55821913e-04, 1.32372938e-04,
#[Out]#                 1.70474435e-04,-7.67997190e-05, 1.25382983e-04,
#[Out]#                 2.19435926e-04, 2.49170873e-04, 1.21254088e-04,
#[Out]#                 1.20941004e-04, 1.47097540e-04, 3.86035645e-05,
#[Out]#                 1.99802642e-04, 4.05528641e-04, 3.70344904e-04,
#[Out]#                 3.66640219e-04, 4.92181629e-04, 6.89130567e-04,
#[Out]#                 7.09335087e-04, 1.03754562e-03, 1.46832969e-03,
#[Out]#                 2.34481157e-03, 2.67038727e-03, 3.38229188e-03,
#[Out]#                 4.08489956e-03, 4.16620681e-03, 3.76376626e-03,
#[Out]#                 3.09925643e-03, 2.83725280e-03, 2.95383856e-03,
#[Out]#                 3.21088848e-03, 3.32037383e-03, 3.72184417e-03,
#[Out]#                 4.10958845e-03, 3.76806292e-03, 3.46868439e-03,
#[Out]#                 2.55515263e-03, 2.06674915e-03, 1.43356563e-03,
#[Out]#                 1.23755168e-03, 9.70762048e-04, 6.38304569e-04,
#[Out]#                 8.37306958e-04, 8.41159315e-04, 6.70382404e-04,
#[Out]#                 5.93354460e-04, 7.73442385e-04, 4.36950970e-04,
#[Out]#                 5.02673618e-04, 6.82109734e-04, 6.52675400e-04,
#[Out]#                 5.76858991e-04, 5.32141305e-04, 6.53596770e-04,
#[Out]#                 9.95970564e-04, 1.42732414e-03, 2.09050998e-03,
#[Out]#                 2.49399478e-03, 2.75292457e-03, 2.88200192e-03,
#[Out]#                 3.05112405e-03, 2.98503973e-03, 2.41113175e-03,
#[Out]#                 1.81707658e-03, 1.32427202e-03, 9.49216250e-04,
#[Out]#                 7.13231857e-04, 5.75712940e-04, 6.31175120e-04,
#[Out]#                 5.58406638e-04, 4.07586514e-04, 4.16627736e-04,
#[Out]#                 2.50438636e-04, 1.62659126e-04, 2.11848310e-04,
#[Out]#                 7.63329226e-05, 1.21856312e-04,-5.82935609e-05,
#[Out]#                -1.22712936e-05, 2.12360756e-05, 2.62116431e-04,
#[Out]#                 2.56484025e-04, 2.07223071e-04, 4.55898378e-04,
#[Out]#                 3.15643993e-04, 2.13611464e-04, 2.27090437e-04,
#[Out]#                 8.71263765e-05, 1.17242409e-04, 3.14424455e-04,
#[Out]#                 2.26660879e-04, 2.66332296e-04, 2.45855190e-04,
#[Out]#                 4.34354151e-04, 4.73582302e-04, 4.96958732e-04,
#[Out]#                 6.30042050e-04, 7.81724812e-04, 1.14549766e-03,
#[Out]#                 1.78976543e-03, 2.13206303e-03, 2.63987412e-03,
#[Out]#                 2.91511696e-03, 3.21121002e-03, 3.15929577e-03,
#[Out]#                 2.56777857e-03, 2.01710546e-03, 1.44553941e-03,
#[Out]#                 1.15906727e-03, 7.57485628e-04, 7.69177277e-04,
#[Out]#                 4.87740181e-04, 2.96747952e-04, 2.19628535e-04,
#[Out]#                 2.78775842e-04, 2.75645289e-04, 3.70938040e-04,
#[Out]#                 4.54537192e-04, 4.19201679e-04, 2.18385147e-04,
#[Out]#                 1.29614564e-04, 4.25238424e-04, 4.14301845e-04,
#[Out]#                 2.83032481e-04, 1.47980361e-04, 1.17112804e-05,
#[Out]#                 1.10899418e-04, 2.96211743e-04, 1.80970819e-04,
#[Out]#                 1.34062109e-04, 1.43008932e-04, 1.65949517e-04,
#[Out]#                 8.37772241e-05, 4.01564030e-05,-7.00436576e-05,
#[Out]#                 8.18039189e-05,-8.87141068e-05,-1.71332722e-04,
#[Out]#                -1.69339419e-05,-1.89532962e-04,-2.08710524e-04,
#[Out]#                -1.56719922e-04,-1.56140068e-05,-1.63593781e-04,
#[Out]#                 5.51197518e-05,-3.80688653e-06,-1.38494424e-05,
#[Out]#                 1.46130405e-04, 7.90678750e-05, 8.52478843e-05,
#[Out]#                 5.56897430e-04, 3.00507003e-04, 2.90008538e-05,
#[Out]#                 1.41403842e-04, 5.17629860e-06, 2.58707849e-04,
#[Out]#                 5.59874868e-04, 5.90444717e-04, 5.23205614e-04,
#[Out]#                 8.39197834e-04, 1.08477077e-03, 1.30642985e-03,
#[Out]#                 1.73263450e-03, 1.89991645e-03, 1.84754061e-03,
#[Out]#                 1.76937180e-03, 1.50248501e-03, 1.05029461e-03,
#[Out]#                 6.00320403e-04, 7.49382365e-04, 4.61272313e-04,
#[Out]#                 2.71904661e-04, 4.19568300e-04, 2.22723567e-04,
#[Out]#                 1.76598733e-05, 2.80421169e-04, 3.95415467e-04,
#[Out]#                 2.04585071e-04, 2.00743641e-04, 2.31598780e-04,
#[Out]#                -1.04791870e-04, 1.03162012e-04,-6.14259577e-08,
#[Out]#                -1.11805653e-04,-2.48954399e-04, 7.31191758e-05,
#[Out]#                 1.30930988e-04,-3.11446747e-05, 1.86156685e-05,
#[Out]#                -5.13434752e-05, 1.24426355e-04, 2.47773889e-04,
#[Out]#                 1.07545020e-04, 1.71397318e-04, 2.04212221e-04,
#[Out]#                 4.21747449e-04, 3.10895732e-04, 2.46975542e-04,
#[Out]#                 4.22891200e-04, 6.27139874e-04, 7.75737164e-04,
#[Out]#                 9.38769546e-04, 6.87207154e-04, 4.52987035e-04,
#[Out]#                 4.56237729e-04, 4.64633340e-04, 3.14328383e-04,
#[Out]#                 1.90254170e-04,-3.97860167e-05,-2.19222886e-04,
#[Out]#                -1.52607827e-04, 9.01465683e-05, 2.89350028e-05,
#[Out]#                -2.08304933e-04,-1.54225854e-04, 6.67789027e-06,
#[Out]#                -1.38826130e-04,-3.03000197e-05,-1.90759180e-04,
#[Out]#                 4.99909802e-05,-3.08652525e-05,-4.48897335e-05,
#[Out]#                -3.36920493e-05, 6.85923151e-05, 1.60259879e-04,
#[Out]#                 1.85202924e-04, 5.70063603e-05,-4.47251587e-05,
#[Out]#                 2.14954023e-04, 4.49081817e-05,-4.76416899e-05,
#[Out]#                -1.19645465e-05, 3.35692544e-04, 1.74265588e-04,
#[Out]#                 4.34257090e-04, 3.95372161e-04, 2.58240354e-04,
#[Out]#                 4.90755658e-04, 6.99470285e-04, 1.23837101e-03,
#[Out]#                 1.47992256e-03, 1.61014567e-03, 1.64462638e-03,
#[Out]#                 1.47521077e-03, 1.42088719e-03, 1.18583592e-03,
#[Out]#                 1.00416970e-03, 1.18364848e-03, 1.08337775e-03,
#[Out]#                 8.81330518e-04, 7.93086772e-04, 8.07144854e-04,
#[Out]#                 5.45976509e-04, 5.83473593e-04, 8.16726591e-04,
#[Out]#                 9.60087113e-04, 1.01650797e-03, 1.08822365e-03,
#[Out]#                 1.01423985e-03, 6.86424028e-04, 3.57120909e-04,
#[Out]#                 2.04664830e-04, 3.64647654e-04, 2.71239463e-04,
#[Out]#                 1.50559072e-06, 1.33575551e-04, 1.51812928e-04,
#[Out]#                 5.14518251e-05, 1.18714561e-04,-1.15870949e-04,
#[Out]#                 5.22774535e-05,-4.23435376e-05,-1.50885229e-04,
#[Out]#                -2.24501418e-04, 4.75610977e-05, 1.22564379e-04,
#[Out]#                 2.38242108e-04, 4.28370695e-04, 4.32706001e-04,
#[Out]#                 2.83839152e-04, 6.01077743e-04, 7.96573469e-04,
#[Out]#                 1.01570413e-03, 9.81742749e-04, 6.88043889e-04,
#[Out]#                 6.48913556e-04, 1.40556163e-04, 1.92903695e-04,
#[Out]#                 1.31114037e-04,-3.82248581e-05,-1.18479227e-04,
#[Out]#                -2.55560990e-05, 1.69330406e-06, 2.51071579e-05,
#[Out]#                -3.19862738e-04,-3.68404108e-05,-1.03954306e-04,
#[Out]#                -9.21254614e-05,-2.70927325e-04,-7.81964482e-05,
#[Out]#                -1.67536753e-04,-1.19310695e-04,-1.09858856e-05,
#[Out]#                 2.48811179e-04, 9.35838470e-05, 5.19757159e-05,
#[Out]#                -5.74863298e-05,-6.42456871e-05,-2.96221024e-05,
#[Out]#                -8.72280507e-05,-2.06058045e-04,-3.14810488e-04,
#[Out]#                -9.09204537e-05, 1.06511681e-04, 1.02920618e-04,
#[Out]#                 1.27778505e-04, 1.78420902e-04, 2.26368837e-04,
#[Out]#                 4.21094068e-04, 4.12208989e-04, 6.30038208e-04,
#[Out]#                 1.19347859e-03, 1.20092207e-03, 1.06896227e-03,
#[Out]#                 1.16729503e-03, 1.10387790e-03, 1.04523997e-03,
#[Out]#                 7.47943646e-04, 8.33257567e-04, 7.62085663e-04,
#[Out]#                 5.20078000e-04, 5.53958933e-04, 7.07598461e-04,
#[Out]#                 9.30643466e-04, 1.00634189e-03, 1.28967804e-03,
#[Out]#                 1.47840544e-03, 1.53257267e-03, 1.54061278e-03,
#[Out]#                 1.33699598e-03, 1.23577868e-03, 8.75584490e-04,
#[Out]#                 4.06782841e-04, 3.77406745e-04, 3.31688090e-04,
#[Out]#                 2.39838642e-04, 2.70603865e-04, 2.67107622e-04,
#[Out]#                 3.46003508e-04, 3.13339871e-04, 3.71151167e-04,
#[Out]#                 2.74533319e-04, 3.15965095e-04, 4.94775304e-04,
#[Out]#                 5.69709344e-04, 4.62741707e-04, 4.78643313e-04,
#[Out]#                 6.02392014e-04, 8.26906122e-04, 5.65609895e-04,
#[Out]#                 4.59492760e-04, 3.14566831e-04, 3.92583839e-04,
#[Out]#                 1.91304003e-04, 2.97684950e-04, 7.96525637e-05,
#[Out]#                -1.48153204e-05, 5.16415785e-05, 1.01257909e-04,
#[Out]#                 1.50112479e-04, 1.16148287e-04, 2.02975294e-04,
#[Out]#                 1.16103780e-04, 1.82191608e-04,-4.81222332e-06,
#[Out]#                 6.91104069e-05,-1.90327395e-04,-1.52903565e-04,
#[Out]#                 4.42005767e-05, 9.69187677e-05, 2.31434184e-04,
#[Out]#                 4.00983496e-04, 2.94295460e-04, 3.44790402e-04,
#[Out]#                 4.10393899e-04, 5.01532457e-04, 5.56372164e-04,
#[Out]#                 3.72272887e-04, 4.87906538e-04, 4.44944424e-04,
#[Out]#                 4.02951962e-04, 2.02591080e-04, 1.49044616e-04,
#[Out]#                 1.27661115e-04, 5.99874911e-05,-1.00019693e-04,
#[Out]#                -2.59100871e-05,-1.94829510e-04,-1.47089086e-04,
#[Out]#                 1.60605661e-04, 8.31881916e-05, 2.29456433e-04,
#[Out]#                 9.73874630e-05, 5.84457121e-05, 1.65178033e-04,
#[Out]#                 2.98045023e-04, 2.15814114e-04, 1.16279371e-04,
#[Out]#                 8.90414158e-05, 2.10845959e-04, 5.83784495e-05,
#[Out]#                 8.17172695e-05, 2.44490628e-04,-7.82338393e-05,
#[Out]#                -1.78350965e-05, 1.56564711e-04,-4.72961146e-05,
#[Out]#                -2.09929582e-04, 2.23488198e-04, 3.24796449e-04,
#[Out]#                 7.43130222e-05,-7.03192563e-05,-1.24180224e-05,
#[Out]#                 1.36093266e-04, 2.41561080e-04, 8.43094313e-05,
#[Out]#                 1.51184053e-04, 2.87703500e-04,-9.45837819e-05,
#[Out]#                -2.17890018e-04,-1.08320103e-06,-1.17823889e-04,
#[Out]#                -9.85356601e-05,-5.24363932e-05,-6.09340605e-05,
#[Out]#                -2.67630230e-06,-4.76217174e-05,-1.21937250e-04,
#[Out]#                -9.07802605e-05,-4.12804457e-06,-2.16574132e-04,
#[Out]#                -2.71698140e-04,-1.31722691e-05,-2.96351373e-05,
#[Out]#                 6.82440532e-06, 4.51813627e-04, 3.54237185e-04,
#[Out]#                 6.34846394e-04, 7.14353519e-04, 1.08816591e-03,
#[Out]#                 7.42393662e-04, 9.17292375e-04, 7.51569110e-04,
#[Out]#                 5.19763213e-04, 2.12915242e-04, 2.40633861e-04,
#[Out]#                 1.08968874e-04, 5.19954228e-05, 2.54020764e-04,
#[Out]#                 5.12195402e-05, 1.82315809e-04, 2.33135375e-04,
#[Out]#                 3.89331981e-04, 3.02001077e-04, 3.05517984e-04,
#[Out]#                 1.96867157e-04, 1.73269756e-04, 2.07212433e-04,
#[Out]#                 1.20777535e-04,-6.27947738e-05,-1.48246778e-04,
#[Out]#                 1.02203849e-04,-4.32623801e-06,-9.35979915e-05,
#[Out]#                -8.61859880e-05, 1.93296451e-04, 2.24187432e-04,
#[Out]#                 1.05054525e-04,-2.13225503e-04,-8.51365985e-05,
#[Out]#                -1.76886897e-04,-3.85383639e-04,-1.93028754e-04,
#[Out]#                -2.01100920e-04,-1.57038274e-04,-3.61566867e-06,
#[Out]#                -1.98310343e-04,-2.36516542e-04,-2.34401014e-04,
#[Out]#                -1.34541930e-04,-3.08209623e-04,-4.32246248e-04,
#[Out]#                -2.23963783e-04,-1.06990992e-05,-2.31081227e-04,
#[Out]#                -2.03293806e-04,-3.27011076e-04,-1.21933794e-04,
#[Out]#                -7.67364545e-05,-1.24834245e-04,-1.08623462e-04,
#[Out]#                -3.49911075e-04,-3.26962356e-04,-5.09796257e-04,
#[Out]#                -3.00224521e-04,-2.46884621e-04,-1.62139346e-04,
#[Out]#                -2.89760064e-04,-3.55782977e-04,-1.28771993e-04,
#[Out]#                -9.94950096e-05,-5.31733276e-05, 9.10480958e-06,
#[Out]#                -1.65692894e-04,-8.76266058e-05, 3.16766236e-05,
#[Out]#                -9.82193815e-05,-1.94886714e-04,-1.91533036e-04,
#[Out]#                -8.73597455e-05,-1.76321755e-05,-1.35849070e-04,
#[Out]#                -2.34548541e-04,-2.59169581e-04,-5.00883580e-05,
#[Out]#                -1.11196081e-04,-2.64765928e-04,-5.40040783e-04,
#[Out]#                -4.14940441e-04,-2.16660323e-04,-1.90058508e-05,
#[Out]#                -1.08763539e-04,-9.97315583e-05,-1.83810567e-04,
#[Out]#                -1.31259905e-04,-2.35278378e-04,-3.05311754e-04,
#[Out]#                -5.74476580e-05,-1.43234269e-04,-1.09065440e-04,
#[Out]#                -1.08422086e-04, 7.55076835e-05, 1.62527387e-04,
#[Out]#                -8.61445096e-06,-3.01045247e-06,-1.81690644e-04,
#[Out]#                -1.58292821e-07,-1.65711208e-05,-2.37481868e-06,
#[Out]#                -9.26057182e-05, 1.07654407e-04,-4.84049406e-06,
#[Out]#                -2.35013445e-04, 1.07149390e-05, 1.00879639e-04,
#[Out]#                 2.12234867e-04, 1.70421423e-04, 1.14056493e-04,
#[Out]#                 1.83024400e-04, 1.65528327e-04,-1.03764134e-04,
#[Out]#                 1.13819893e-04,-1.68229177e-04,-2.56530620e-04,
#[Out]#                -6.84598403e-04,-1.07968051e-03,-2.57667532e-04,
#[Out]#                 1.03715875e-05, 5.72126846e-05, 3.25452682e-04,
#[Out]#                 1.45909158e-04, 1.05716048e-04, 3.01074964e-04,
#[Out]#                 2.33328319e-04, 2.45886185e-04, 1.95644432e-04,
#[Out]#                 1.63090415e-04, 4.86211298e-04, 5.54586179e-04,
#[Out]#                 7.35139765e-04, 3.62230290e-04, 4.35043155e-04,
#[Out]#                 4.90662875e-04, 2.59351102e-04, 2.12888306e-04,
#[Out]#                 1.16278534e-04, 2.77724292e-04, 6.08138507e-05,
#[Out]#                 4.09789900e-05, 3.00492474e-07,-5.95822567e-05,
#[Out]#                -3.62104533e-04,-5.22988266e-04,-3.11623939e-04,
#[Out]#                 2.44722614e-05, 1.30757064e-04,-1.03670078e-04,
#[Out]#                 6.95310373e-05, 3.33714212e-04, 6.24288441e-05,
#[Out]#                -6.78124525e-06,-7.36912014e-04,-2.10027653e-03,
#[Out]#                -1.84046489e-03,-1.33847829e-03,-1.37599220e-03,
#[Out]#                -1.57984055e-03,-1.54615426e-03,-1.69899524e-03,
#[Out]#                -1.31153886e-03,-1.34037330e-03,-1.40292430e-03,
#[Out]#                -1.38526654e-03,-1.71004550e-03,-1.35570997e-03,
#[Out]#                -7.20232725e-04,-3.88476212e-04,-2.38356079e-04,
#[Out]#                 1.26905259e-04,-9.88837419e-05,-6.80096637e-05,
#[Out]#                -1.25864084e-04, 2.12500876e-04, 1.74756904e-04,
#[Out]#                 1.62106531e-04, 9.82640995e-05, 1.02647195e-04,
#[Out]#                 4.36998344e-05,-1.35919836e-04,-3.78334720e-04,
#[Out]#                -1.28290022e-03,-1.83164829e-03,-2.10026349e-03,
#[Out]#                -1.33115123e-03,-4.81046969e-04,-2.44524708e-04,
#[Out]#                -3.68830399e-04,-3.57516663e-04,-5.36548265e-04,
#[Out]#                -4.69899562e-04,-5.67621668e-04,-8.50739656e-04,
#[Out]#                -1.30707864e-03,-1.72491407e-03,-1.39913289e-03,
#[Out]#                -1.31526287e-03] Jy / beam>
djtok = davgspec.beams.common_beam().jtok(davgspec.spectral_axis).mean()
djtok = davgspec.beam.jtok(davgspec.spectral_axis).mean()
djtok = davgspec.beam.jtok(davgspec.with_spectral_unit(u.GHz).spectral_axis).mean()
dsp = pyspeckit.Spectrum(data=(davgspec - np.median(davgspec.quantity))*djtok/davgspec.unit, xarr=davgspec.spectral_axis)
dsp.plotter()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 16, 1], use_lmfit=True)
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.02):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(xarr, vcen, width, tex, column,
                                       freqs, aij, deg, EU, partfunc)*fillingfactor
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d0bc5fca0>]
fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.012], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
dcube = SpectralCube.read('/orange/adamginsburg/brick_alma_linesurvey/danwalker/Brick_CH3CN.image.fits')
reg = regions.Regions.read('../centralcoreellipse.reg')
dscube = dcube.subcube_from_regions(reg)
dscube
#[Out]# SpectralCube with shape=(700, 117, 114) and unit=Jy / beam:
#[Out]#  n_x:    114  type_x: RA---SIN  unit_x: deg    range:   266.543732 deg:  266.544806 deg
#[Out]#  n_y:    117  type_y: DEC--SIN  unit_y: deg    range:   -28.705359 deg:  -28.704392 deg
#[Out]#  n_s:    700  type_s: VRAD      unit_s: m / s  range:   -63094.684 m / s:  400457.822 m / s
m0 = dscube.with_spectral_unit(u.GHz).spectral_slab(220.71*u.GHz, 220.72*u.GHz).moment0() # slab of comps 0 and 1
m0.quicklook()
davgspec = (dscube * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
djtok = davgspec.beam.jtok(davgspec.with_spectral_unit(u.GHz).spectral_axis).mean()
dsp = pyspeckit.Spectrum(data=(davgspec - np.median(davgspec.quantity))*djtok/davgspec.unit, xarr=davgspec.spectral_axis)
dsp.plotter()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 16, 1], use_lmfit=True)
dsp = pyspeckit.Spectrum(data=(davgspec - np.median(davgspec.quantity))*djtok/davgspec.unit, xarr=davgspec.with_spectral_unit(u.GHz).spectral_axis)
dsp.plotter()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e6, 1], use_lmfit=True)
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = 220.3*u.GHz, fmax = 220.9*u.GHz)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.02):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(xarr, vcen, width, tex, column,
                                       freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = 220.3*u.GHz, fmax = 220.9*u.GHz)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=1):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(xarr, vcen, width, tex, column,
                                       freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
m0.argmax()
#[Out]# 0
np.nanargmax(m0)
#[Out]# 5528
np.unravel_index(np.nanargmax(m0), m0.shape)
#[Out]# (48, 56)
#davgspec = (dscube * m0.value).sum(axis=(1,2)) / np.nansum(m0.value)
davgspec = dscube[:,48,56]
djtok = davgspec.beam.jtok(davgspec.with_spectral_unit(u.GHz).spectral_axis).mean()
dsp = pyspeckit.Spectrum(data=(davgspec - np.median(davgspec.quantity))*djtok/davgspec.unit, xarr=davgspec.with_spectral_unit(u.GHz).spectral_axis)
dsp.plotter()
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = 220.3*u.GHz, fmax = 220.9*u.GHz)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=1):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(xarr, vcen, width, tex, column,
                                       freqs, aij, deg, EU, partfunc)*fillingfactor
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

dsp.specfit.Registry.multifitters['CH3CN'] = fitter
dsp.specfit.Registry.npars['CH3CN'] = 5

dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 195.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
dsp.plotter()
dsp.specfit.plot_fit()
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 195.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 105.1, 1e16, 1], use_lmfit=True)
dsp.plotter()
dsp.specfit.plot_fit()
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
get_ipython().run_line_magic('pinfo2', 'dsp.specfit.plot_model')
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 3.6, 167, 8.9e16])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 3.6, 167, 8.9e16, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 3.6, 167, 8.9e15, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 3.6, 167, 8.9e14, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 3.6, 167, 2e15, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 2.6, 167, 2e15, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 2.6, 167, 1e15, 1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 2.6, 167, 8.6e16, 0.2])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 2.6, 167, 8.6e16, 0.1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
dsp.plotter()
dsp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.1])
dsp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
geg
deg
#[Out]# array([ 50, 100,  50,  50, 100,  50,  50, 100,  50,  50,  50, 100,  50])
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.1])
sp.plotter.axis.set_ylim(-5,30)
#[Out]# (-5.0, 30.0)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.1])
sp.plotter.axis.set_ylim(-1,3)
#[Out]# (-1.0, 3.0)
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.model import SpectralModel
sp = pyspeckit.Spectrum(data=(avgspec - np.median(avgspec.quantity))*jtok/avgspec.unit, xarr=avgspec.spectral_axis)
sp.plotter()
v_cen = 40*u.km/u.s
v_disp = 2*u.km/u.s
temp = 167*u.K
N_tot = 1e16*u.cm**-2
species = 'CH3CN'
freqs, aij, deg, EU, partfunc = get_molecular_parameters(species, fmin = fmin, fmax = fmax)
mod = lte_molecule.generate_model(sp.xarr, v_cen, v_disp, temp, N_tot,
                                  freqs, aij, deg, EU, partfunc)
def modfunc(xarr, vcen, width, tex, column, fillingfactor=0.02):
    if column < 100:
        column = 10**column
    return lte_molecule.generate_model(xarr, vcen, width, tex, column,
                                       freqs, aij, deg, EU, partfunc)*fillingfactor
sp.xarr.convert_to_unit(u.km/u.s, refX=freqs[-1], velocity_convention='radio')
mod = modfunc(sp.xarr, v_cen.value, v_disp.value, temp.value, N_tot.value)
sp.plotter()
sp.plotter.axis.plot(sp.xarr, mod)
#[Out]# [<matplotlib.lines.Line2D at 0x2b3d12142c40>]
fitter = SpectralModel(modfunc, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 4

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 1e15], debug=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
fitter = SpectralModel(modfunc, 5,
            parnames=['shift','width','tex','column','ff'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,True)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,1)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','ff'),
            centroid_par='shift',
            )
fitter.__name__ = "CH3CN"

sp.specfit.Registry.multifitters['CH3CN'] = fitter
sp.specfit.Registry.npars['CH3CN'] = 5

sp.specfit(fittype='CH3CN', guesses=[40.1, 1.5, 155.1, 15, 0.012], use_lmfit=True)
sp.plotter()
sp.specfit.plot_fit()
sp.specfit.plotresiduals()
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.01])
sp.plotter.axis.set_ylim(-1,3)
#[Out]# (-1.0, 3.0)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.001])
sp.plotter.axis.set_ylim(-1,3)
#[Out]# (-1.0, 3.0)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.001])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e16, 0.005])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e15, 0.05])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 2.6, 167, 8.9e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 167, 8.9e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 167, 3e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 167, 6e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 107, 6e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 57, 6e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 57, 3e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 57, 1e15, 0.02])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 57, 1e15, 0.05])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
sp.plotter()
sp.specfit.plot_model([41.05, 1.6, 57, 1e15, 0.04])
sp.plotter.axis.set_ylim(-1,1.5)
#[Out]# (-1.0, 1.5)
