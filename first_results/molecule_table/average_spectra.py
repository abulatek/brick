### This script was used to extract spectra from continuum-subtracted smoothed cubes; 
### see average_spectra.ipynb for a script that was used to extract spectra from NON-contsub cubes.
# sbatch --job-name=extract_contsub_spectra --output=extract_contsub_spectra-%j.log --account=astronomy-dept --qos=astronomy-dept-b --ntasks=8 --mail-type=ALL --mail-user=abulatek@ufl.edu --nodes=1 --mem=24gb --time=96:00:00 --constraint=el8 --wrap "/blue/adamginsburg/abulatek/miniconda3/bin/python3 average_spectra.py"

# Package imports
from astropy import coordinates, units as u
import warnings
warnings.filterwarnings('ignore')
from astropy import log
log.setLevel('ERROR')
from spectral_cube import SpectralCube
import numpy as np

# My functions
from mol_model import fetch_cubes, model_and_plot, list_mol_tags, get_cubes_from_mask

results = '/blue/adamginsburg/abulatek/brick/symlinks/contsub_smoothed_cubes_K/'
avg_spectra_core_loc = '/blue/adamginsburg/abulatek/brick/symlinks/contsub_smoothed_cubes_K/avg_spectra_core/'
avg_spectra_frown_loc = '/blue/adamginsburg/abulatek/brick/symlinks/contsub_smoothed_cubes_K/avg_spectra_frown/'

# Get list of all cubes
import glob
cubefns = glob.glob(f"{results}/source_ab_*.image.commonbeam.contsub.K.fits")

# Get all cubes in order
cubes = []
for fn in cubefns:
    molcube = SpectralCube.read(fn, format='fits')
    cubes.append(molcube)
# Reorder cubes from lowest to highest frequency
lowestfreqs = [np.min(molcube.spectral_axis).value for molcube in cubes]
correct_order = np.argsort(lowestfreqs)
cubes_inorder = [cubes[i] for i in correct_order]

# Set coordinate of central source
crd = coordinates.SkyCoord("17:46:10.6339473267 -28:42:17.9807702398", frame='icrs', unit=(u.h, u.deg))

# Get all spectra in order for central source coordinate
spectra = []
for molcube in cubes_inorder:
    # Extract spectrum from provided coordinate
    x,y = map(int, molcube.wcs.celestial.world_to_pixel(crd))
    data_sp = molcube[:, y, x]
    # Continuum-subtract the spectrum
    data_sp_contsub = data_sp - np.median(data_sp)
    spectra.append(data_sp_contsub)
    
# Save spectra for central source coordinate
for spectrum in spectra:
    fn = f"spectrum_{np.min(spectrum.spectral_axis).value}_{np.max(spectrum.spectral_axis).value}.fits"
    spectrum.write(avg_spectra_core_loc+fn, format='fits')
    
# Get subcubes based on mask
cubes_inorder_masked = get_cubes_from_mask("diffuse_regions.fits", 1, cubes_inorder, plot_region=False)

# Get all spectra in order for frown region
spectra_reg = []
for molcube in cubes_inorder_masked:
    # Extract spectrum from the subcube made from the mask
    data_sp_reg = molcube.mean(axis=(1,2), how='slice')
    # Continuum-subtract the spectrum
    data_sp_reg_contsub = data_sp_reg - np.median(data_sp_reg)
    spectra_reg.append(data_sp_reg_contsub)

# Save spectra for frown region
for spectrum in spectra_reg:
    fn = f"spectrum_{np.min(spectrum.spectral_axis).value}_{np.max(spectrum.spectral_axis).value}.fits"
    spectrum.write(avg_spectra_frown_loc+fn, format='fits')