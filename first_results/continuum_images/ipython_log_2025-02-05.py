########################################################
# Started Logging At: 2025-02-05 14:07:07
########################################################
########################################################
# # Started Logging At: 2025-02-05 14:07:08
########################################################
########################################################
# Started Logging At: 2025-02-05 14:07:27
########################################################
########################################################
# # Started Logging At: 2025-02-05 14:07:27
########################################################
########################################################
# Started Logging At: 2025-02-05 14:07:38
########################################################
########################################################
# # Started Logging At: 2025-02-05 14:07:38
########################################################
from spectral_cube import SpectralCube
import spectral_cube
print(spectral_cube.__version__)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
from astropy.io import fits
from reproject import reproject_interp
from astropy import wcs, units as u
results = '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/'
# Get list of all images
import glob
imagefns = glob.glob(f"{results}/brickmaser_cont_*_mtmfs_incl_adjparam*.image.tt0")
# Remove the whole-bandwidth versions of data_31 and data_41
imagefns.remove('/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_31_mtmfs_incl_adjparam.image.tt0')
imagefns.remove('/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_41_mtmfs_incl_adjparam.image.tt0')
imagefns
#[Out]# ['/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_41_mtmfs_incl_adjparam_split_first.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_31_mtmfs_incl_adjparam_split_second.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_55_mtmfs_incl_adjparam.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_59_mtmfs_incl_adjparam.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_51_mtmfs_incl_adjparam.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_41_mtmfs_incl_adjparam_split_second.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_5d_mtmfs_incl_adjparam.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_31_mtmfs_incl_adjparam_split_first.image.tt0',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/continuum_images/brickmaser_cont_data_61_mtmfs_incl_adjparam.image.tt0']
# Import the continuum images
contimages = []
for imagefn in imagefns:
    image = SpectralCube.read(imagefn, format='casa_image')
    # image = image.to(u.K) # We'd rather have Jy/beam because we're dealing with point sources
    contimages.append(image)
# Reorder continuum images from lowest to highest frequency
lowestfreqs = [np.min(contimage.spectral_axis).value for contimage in contimages]
correct_order = np.argsort(lowestfreqs)
contimages_ordered = [contimages[i] for i in correct_order]
# Calculate common beam
import radio_beam
beams = [contimage.beam for contimage in contimages_ordered]
beams = radio_beam.Beams(beams=beams)
common_beam = beams.common_beam()
common_beam
# Convert all images to a common beamsize
contimages_ordered_cb = [contimage.convolve_to(common_beam) for contimage in contimages_ordered]
# # Create a mask for the maser core (only have to run this once)
# mask = contimages_ordered[5][0] > 2.*u.mJy/u.beam
# mask_path = "/blue/adamginsburg/abulatek/brick/first_results/continuum_images/maser_core.fits"
# fits.PrimaryHDU(data=mask.astype('int'), header=contimages_ordered[5][0].wcs.celestial.to_header()).writeto(mask_path, overwrite=True)
# Set core region
core_reg = fits.open('maser_core.fits')
core_reg_mask = core_reg[0].data == 1
test_map = contimages_ordered_cb[0].with_mask(core_reg_mask)
plt.imshow(test_map[0].value, origin="lower")
plt.show()
# Calculate the peak intensities in each of the continuum images (with central core mask on)
freqs = []
peaks = []
for contimage in contimages_ordered_cb:
    freq = contimage.spectral_axis[0].to(u.GHz) # This is the center freq of the ms
    # peak = np.nanmax(contimage[0].value)*contimage.unit # Peak intensity of the whole image
    peak = np.nanmax(contimage.with_mask(core_reg_mask)[0].value)*contimage[0].unit # Peak intensity of core region
    freqs.append(freq)
    peaks.append(peak)
freqs_vals = np.array([freq.value for freq in freqs])
peaks_vals = np.array([peak.value for peak in peaks])
# for contimage in contimages_ordered_cb:
#     immin = np.nanmin(contimage[0].value)
#     immax = np.nanmax(contimage[0].value)
#     ax = plt.subplot(111, projection = contimage[0].wcs)
#     im = ax.imshow(contimage[0].value, origin = 'lower', cmap='viridis', norm='linear', vmax=0.11*immax, vmin=immin)
#     cbar = plt.colorbar(im, ax=ax)
#     cbar.set_label(f"Brightness temperature [{contimage[0].unit.to_string(format = 'latex_inline')}]")
#     ax.set_ylabel('Declination')
#     ax.set_xlabel('Right ascension') 
#     # Add core contour
#     mask = contimages_ordered[5][0] > 2.*u.mJy/u.beam
#     ax.contour(mask, levels = [0, 1], linewidths=0.75, colors = ['r'])
#     # Put beam on each image
#     pixscale = np.abs((wcs.utils.proj_plane_pixel_area(contimage[0].wcs)**0.5*u.deg).to(u.arcsec))
#     bm = contimage[0].beam
#     bmell = bm.ellipse_to_plot(15, 15, pixscale)
#     bmell.set_facecolor('none')
#     bmell.set_edgecolor('k')
#     ax.add_artist(bmell)
#     # Add center frequency to plot
#     freq = str(round(contimage.spectral_extrema[0].to(u.GHz).value, 2))+' GHz'
#     plt.text(0.02, 0.95, freq, transform=ax.transAxes)
#     # Compute the angle corresponding to some distance in parsecs at the distance of the galactic center
#     from astropy.visualization.wcsaxes import add_scalebar
#     gc_distance = 8.0*u.kpc
#     scalebar_length = 0.25*u.pc
#     scalebar_angle = (scalebar_length / gc_distance).to(u.deg, equivalencies=u.dimensionless_angles())
#     # Add a scale bar
#     add_scalebar(ax, scalebar_angle, label=f"{scalebar_length.value} {scalebar_length.unit}", color="black", frame=True)
#     plt.show()
from scipy.optimize import curve_fit
def linear(x, m, b):
    return m*x + b
# Do a linear fit on the core SED
popt, pcov = curve_fit(linear, freqs_vals, peaks_vals)
slope, intercept = popt[0], popt[1]
# Import the continuum images
contimages = []
for imagefn in imagefns:
    image = SpectralCube.read(imagefn, format='casa_image')
    # image = image.to(u.K) # We'd rather have Jy/beam because we're dealing with point sources
    contimages.append(image)
# Reorder continuum images from lowest to highest frequency
lowestfreqs = [np.min(contimage.spectral_axis).value for contimage in contimages]
correct_order = np.argsort(lowestfreqs)
contimages_ordered = [contimages[i] for i in correct_order]
# Calculate common beam
import radio_beam
beams = [contimage.beam for contimage in contimages_ordered]
beams = radio_beam.Beams(beams=beams)
common_beam = beams.common_beam()
common_beam
# Convert all images to a common beamsize
contimages_ordered_cb = [contimage.convolve_to(common_beam) for contimage in contimages_ordered]
for imagefn in imagefns:
    image = SpectralCube.read(imagefn, format='casa_image')
    conv = image.convolve_to(commonbeam)
    conv.write(imagefn+".commonbeam.fits")
# Import the continuum images
contimages = []
for imagefn in imagefns:
    image = SpectralCube.read(imagefn, format='casa_image')
    # image = image.to(u.K) # We'd rather have Jy/beam because we're dealing with point sources
    contimages.append(image)
# Reorder continuum images from lowest to highest frequency
lowestfreqs = [np.min(contimage.spectral_axis).value for contimage in contimages]
correct_order = np.argsort(lowestfreqs)
contimages_ordered = [contimages[i] for i in correct_order]
# Calculate common beam
import radio_beam
beams = [contimage.beam for contimage in contimages_ordered]
beams = radio_beam.Beams(beams=beams)
common_beam = beams.common_beam()
common_beam
# Convert all images to a common beamsize
contimages_ordered_cb = [contimage.convolve_to(common_beam) for contimage in contimages_ordered]
for imagefn in imagefns:
    image = SpectralCube.read(imagefn, format='casa_image')
    conv = image.convolve_to(common_beam)
    conv.write(imagefn+".commonbeam.fits")
