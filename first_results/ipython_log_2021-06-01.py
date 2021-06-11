########################################################
# Started Logging At: 2021-06-01 15:59:06
########################################################
########################################################
# # Started Logging At: 2021-06-01 16:01:31
########################################################
########################################################
# Started Logging At: 2021-06-01 16:05:24
########################################################

########################################################
# # Started Logging At: 2021-06-01 16:05:44
########################################################
from spectral_cube import SpectralCube
# Import cubes
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
fn_87_spw25 = results+'source_ab_87_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
fn_102_spw106 = results+'source_ab_102_spw106_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
fn_95_spw25 = results+'source_ab_95_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image'
cube_87_spw25 = SpectralCube.read(fn_87_spw25, format='casa_image', use_dask=True) 
cube_102_spw106 = SpectralCube.read(fn_102_spw106, format='casa_image', use_dask=True) 
cube_95_spw25 = SpectralCube.read(fn_95_spw25, format='casa_image', use_dask=True)
cube_85_spw25
cube_86_spw25
cube_87_spw25
tbl = get_line(cube_87_spw25, ' H13CN ', 86.3*u.GHz, 86.36*u.GHz, 9.5)
tbl.show_in_notebook()
def get_line(cube, chemical_name, freq_lo, freq_hi, vel):
        tbl = Splatalogue.query_lines(freq_lo-0.1*u.GHz, freq_hi+0.1*u.GHz, 
                                          chemical_name=chemical_name,
                                                                            energy_max=500, # more lines w/ max energy > 140
                                                                                                              energy_type='eu_k',
                                                                                                                                                line_lists=['JPL'],
                                                                                                                                                                                  show_upper_degeneracy=True, 
                                                                                                                                                                                                                    show_qn_code=True)
                                                                                                                                                                                                                        line_freqs = tbl['Meas Freq-GHz(rest frame,redshifted)'].data
shifted_line_freqs = line_freqs-((vel/299792)*line_freqs) # Shift by velocity
tbl['Shifted Freq-GHz'] = shifted_line_freqs
return tbl
def get_subcube(cube, center_freq, slab_width):
        print(center_freq)
            subcube = cube.spectral_slab(center_freq - 0.1*u.GHz, center_freq + 0.1*u.GHz).to(u.K)
subcube_v = subcube.with_spectral_unit(u.km/u.s, 
                                       rest_value = center_freq,
                                                                              velocity_convention = 'radio').spectral_slab(-slab_width,
                                                                                                                                                                  slab_width)
print(subcube_v)
return subcube_v
def get_moments(subcube_v):
        moment_0 = subcube_v.moment0()
            moment_1 = subcube_v.moment1()
moment_2 = subcube_v.moment2()
sigma_map = subcube_v.linewidth_sigma()  
fwhm_map = subcube_v.linewidth_fwhm()  
return moment_0, moment_1, moment_2, sigma_map, fwhm_map
get_ipython().run_line_magic('cpaste', '')
def get_line(cube, chemical_name, freq_lo, freq_hi, vel):
    tbl = Splatalogue.query_lines(freq_lo-0.1*u.GHz, freq_hi+0.1*u.GHz, 
                                  chemical_name=chemical_name,
                                  energy_max=500, # more lines w/ max energy > 140
                                  energy_type='eu_k',
                                  line_lists=['JPL'],
                                  show_upper_degeneracy=True, 
                                  show_qn_code=True)
    line_freqs = tbl['Meas Freq-GHz(rest frame,redshifted)'].data
    shifted_line_freqs = line_freqs-((vel/299792)*line_freqs) # Shift by velocity
    tbl['Shifted Freq-GHz'] = shifted_line_freqs
    return tbl

def get_subcube(cube, center_freq, slab_width):
    print(center_freq)
    subcube = cube.spectral_slab(center_freq - 0.1*u.GHz, center_freq + 0.1*u.GHz).to(u.K)
    subcube_v = subcube.with_spectral_unit(u.km/u.s, 
                                           rest_value = center_freq,
                                           velocity_convention = 'radio').spectral_slab(-slab_width,
                                                                                        slab_width)
    print(subcube_v)
    return subcube_v

def get_moments(subcube_v):
    moment_0 = subcube_v.moment0()
    moment_1 = subcube_v.moment1()
    moment_2 = subcube_v.moment2()
    sigma_map = subcube_v.linewidth_sigma()  
    fwhm_map = subcube_v.linewidth_fwhm()  
    return moment_0, moment_1, moment_2, sigma_map, fwhm_map
tbl = get_line(cube_87_spw25, ' H13CN ', 86.3*u.GHz, 86.36*u.GHz, 9.5)
tbl.show_in_notebook()
from astropy import units as u
tbl = get_line(cube_87_spw25, ' H13CN ', 86.3*u.GHz, 86.36*u.GHz, 9.5)
tbl.show_in_notebook()
from astroquery.splatalogue import Splatalogue
tbl.show_in_notebook()
tbl = get_line(cube_87_spw25, ' H13CN ', 86.3*u.GHz, 86.36*u.GHz, 9.5)
tbl
center_freq = tbl['Shifted Freq-GHz'][-1]*u.GHz
cube_5 = get_subcube(cube_87_spw25, center_freq, 5.*u.km/u.s)
cube_30 = get_subcube(cube_87_spw25, center_freq, 30.*u.km/u.s)
cube_100 = get_subcube(cube_87_spw25, center_freq, 100.*u.km/u.s)
cube_30
center_freq = 102.924253*u.GHz
cube = get_subcube(cube_102_spw106, center_freq, 15.*u.km/u.s)
cube_sclip = cube.sigma_clip_spectrally(3) # Clip values above 3-sigma
mad_std_spectrum_sclip = cube_sclip.mad_std(axis=(1, 2))
mad_std_spectrum_sclip
mad_std_map_sclip = cube_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension
# Make a low and high mask
low_snr_mask = cube_30 > 3 * mad_std_map_sclip
high_snr_mask = cube_30 > 6 * mad_std_map_sclip
high_snr_mask
high_snr_mask.include()
high_snr_mask.include().compute()
# From the labels, count the number of pixels within each label.
num_pixels_in_high_snr_mask = ndmeasure.sum_labels(high_snr_mask.include(),
                                                   label_image=high_snr_mask_labels,
                                                                                                      index=range(1, num_labels.compute() + 1))
from dask_image import ndmeasure
get_ipython().run_line_magic('pip', 'install dask_image')
import dask_image
from dask_image import ndmeasure
structure = np.ones((3, 3, 3), dtype=bool)
low_snr_mask_labels, num_labels = ndmeasure.label(low_snr_mask.include(),
                                                  structure=structure)
print(f"Initial number of regions found: {num_labels.compute()}")
low_snr_mask
cube_30
# From the labels, count the number of pixels within each label.
num_pixels_in_high_snr_mask = ndmeasure.sum_labels(high_snr_mask.include(),
                                                   label_image=high_snr_mask_labels,
                                                                                                      index=range(1, num_labels.compute() + 1))
high_snr_mask_labels, num_labels_high = ndmeasure.label(high_snr_mask.include(),
                                                  structure=structure)
num_labels_high
num_labels_high.compute()
num_pixels_in_high_snr_mask = ndmeasure.sum_labels(high_snr_mask.include(),
                                                   label_image=high_snr_mask_labels,
                                                                                                      index=range(1, num_labels.compute() + 1))
num_pixels_in_high_snr_mask
num_pixels_in_high_snr_mask.compute()
hnm = high_snr_mask.include().compute()
hnm
high_snr_mask_labels
high_snr_mask_labels.compute
from dask.diagnostics import ProgressBar
pbar = ProgressBar()
pbar.register()
high_snr_mask_labels.compute()
num_pixels_in_high_snr_mask = ndmeasure.sum_labels(high_snr_mask.include(),
                                                   label_image=high_snr_mask_labels,
                                                                                                      index=range(1, num_labels.compute() + 1))
num_pixels_in_high_snr_mask
num_pixels_in_high_snr_mask.compute()
# Repeat for the high signal mask.
num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask.include(),
                                                  label_image=low_snr_mask_labels,
                                                                                                    index=range(1, num_labels.compute() + 1))
num_pixels_in_low_snr_mask
num_labels_high.compute()
num_pixels_in_low_snr_mask.compute()
low_snr_mask_labels
lsnl = low_snr_mask_labels.compute()
import time
t0 = time.time(); [(lsnl == label).sum() for label in range(1,6570)]; print(time.time() - t0)
import scipy.ndimage
dir(scipy.ndimage)
scipy.ndimage.sum_labels(lsnl, label_image=low_snr_mask_labels, index=range(1, low_snr_mask_labels.max()+1))
get_ipython().run_line_magic('pinfo', 'scipy.ndimage.sum')
t0=time.time(); scipy.ndimage.sum(lsnl, label_image=low_snr_mask_labels, index=range(1, low_snr_mask_labels.max()+1)); print(time.time() - t0)
lsnl
low_snr_mask_labels
lsnrml = low_snr_mask_labels.compute()
t0=time.time(); scipy.ndimage.sum(lsnl, label_image=lsnrml, index=range(1, lsnrml.max()+1)); print(time.time() - t0)
t0=time.time(); scipy.ndimage.sum(lsnl, image=lsnrml, index=range(1, lsnrml.max()+1)); print(time.time() - t0)
get_ipython().run_line_magic('pinfo', 'scipy.ndimage.sum')
t0=time.time(); scipy.ndimage.sum(lsnl, labels=lsnrml, index=range(1, lsnrml.max()+1)); print(time.time() - t0)
44.14/0.14
get_ipython().run_line_magic('pinfo2', 'ndmeasure.sum_labels')
get_ipython().run_line_magic('pinfo2', 'scipy.ndimage.sum')
get_ipython().run_line_magic('pinfo', 'ndmeasure.sum')
get_ipython().run_line_magic('pinfo2', 'ndmeasure.sum')
get_ipython().run_line_magic('pinfo2', 'ndmeasure.labeled_comprehension')
get_ipython().run_line_magic('pinfo2', 'ndmeasure.sum_labels')
get_ipython().run_line_magic('pinfo', 'ndmeasure.sum_labels')
get_ipython().run_line_magic('pinfo2', 'ndmeasure.sum_labels')
ndmeasure._utils
get_ipython().run_line_magic('pinfo', 'ndmeasure._utils._histogram')
get_ipython().run_line_magic('pinfo2', 'ndmeasure._utils._norm_input_labels_index')
get_ipython().run_line_magic('pinfo2', 'ndmeasure.labeled_comprehension')
get_ipython().run_line_magic('pinfo2', 'ndmeasure._utils._labeled_comprehension_func')
get_ipython().run_line_magic('pinfo2', 'ndmeasure._utils._labeled_comprehension_delayed')
get_ipython().run_line_magic('pinfo2', 'ndmeasure._utils._labeled_comprehension_func')
get_ipython().run_line_magic('pinfo2', 'ndmeasure._utils._labeled_comprehension_delayed')
get_ipython().run_line_magic('history', '')
lsnrml
low_snr_mask
low_snr_mask_labels
low_snr_mask_labels.rechunk(19,256,256)
low_snr_mask_labels.rechunk((19,256,256))
get_ipython().run_line_magic('history', '')
#num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask.include(), label_image=low_snr_mask_labels, index=range(1, num_labels.compute() + 1))
low_snr_mask.include()
num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask.include().rechunk((19,512,512)), label_image=low_snr_mask_labels.rechunk((19,512,512)), index=range(1, num_labels.compute() + 1))
t0=time.time(); num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask.include().rechunk((19,512,512)), label_image=low_snr_mask_labels.rechunk((19,512,512)), index=range(1, num_labels.compute() + 1)); print(time.time() - t0)
t0=time.time(); lsnrml = low_snr_mask_labels.compute(); lsnl = low_snr_mask_labels.compute(); rslt = scipy.ndimage.sum(lsnl, label_image=low_snr_mask_labels, index=range(1, low_snr_mask_labels.max()+1)); print(time.time() - t0)
t0=time.time(); lsnrml = low_snr_mask_labels.compute(); lsnl = low_snr_mask_labels.compute(); rslt = scipy.ndimage.sum(lsnl, label_image=lsnrml, index=range(1, lsnrml.max()+1)); print(time.time() - t0)
t0=time.time(); lsnrml = low_snr_mask_labels.compute(); lsnl = low_snr_mask_labels.compute(); rslt = scipy.ndimage.sum(lsnl, labels=lsnrml, index=range(1, lsnrml.max()+1)); print(time.time() - t0)
t0=time.time(); num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask.include().rechunk((19,512,512)), label_image=low_snr_mask_labels.rechunk((19,512,512)), index=range(1, num_labels.compute() + 1)); print(time.time() - t0)
blah = np.random.randn(19,512,512)
blah.shape
msk = blah > 1
scipy.ndimage.label(msk)
lab, ct = scipy.ndimage.label(msk)
ct
msk = blah > 2
lab, ct = scipy.ndimage.label(msk)
ct
msk = blah > 3
lab, ct = scipy.ndimage.label(msk)
ct
scipy.ndimage.sum(msk, labels=lab, index=range(1, ct+1))
get_ipython().run_line_magic('timeit', 'scipy.ndimage.sum(msk, labels=lab, index=range(1, ct+1))')
get_ipython().run_line_magic('history', '')
ndmeasure.sum_labels(msk, label_image=lab, index=range(1, ct+1))
rslt = ndmeasure.sum_labels(msk, label_image=lab, index=range(1, ct+1))
rslt.compute()
rslt.rechunk(1000).compute()
22.8 / 0.117
get_ipython().run_line_magic('timeit', 'ndmeasure.sum_labels(msk, label_image=lab, index=range(1, ct+1))')
