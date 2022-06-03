########################################################
# Started Logging At: 2022-05-18 15:18:16
########################################################
########################################################
# # Started Logging At: 2022-05-18 15:18:17
########################################################
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'

# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
cubefns = glob.glob(f"{results}/*.image.pbcor")
get_ipython().run_line_magic('ls', '$results/*pbcor* -d')
from astropy import coordinates, units as u
from astropy import wcs
from astropy.io import fits
crd = coordinates.SkyCoord("17:46:10.80 -28:42:13.3", frame='fk5', unit=(u.h, u.deg))
import warnings
warnings.filterwarnings('ignore')
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn, format='casa_image')
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
        print(fn)
        print(ch3cntbl)
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
        ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.pbcor.fits', overwrite=True)
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
import numpy as np
# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn, format='casa_image')
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
        print(fn)
        print(ch3cntbl)
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
        ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.pbcor.fits', overwrite=True)
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(59765, 336, 336) and unit=Jy / beam and chunk size (934, 24, 336):
#[Out]#  n_x:    336  type_x: RA---SIN  unit_x: deg    range:   266.528944 deg:  266.559715 deg
#[Out]#  n_y:    336  type_y: DEC--SIN  unit_y: deg    range:   -28.718463 deg:  -28.691477 deg
#[Out]#  n_s:  59765  type_s: FREQ      unit_s: Hz     range: 125087368295.955 Hz:154269854424.339 Hz
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn, format='casa_image')
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
        print(fn)
        print(ch3cntbl)
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        for rf, qn in zip(restfrq, qns):
            print(fn, qn, rf)
            ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qn}.pbcor.fits', overwrite=True)
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(59765, 336, 336) and unit=Jy / beam and chunk size (934, 24, 336):
#[Out]#  n_x:    336  type_x: RA---SIN  unit_x: deg    range:   266.528944 deg:  266.559715 deg
#[Out]#  n_y:    336  type_y: DEC--SIN  unit_y: deg    range:   -28.718463 deg:  -28.691477 deg
#[Out]#  n_s:  59765  type_s: FREQ      unit_s: Hz     range: 125087368295.955 Hz:154269854424.339 Hz
tempmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
colmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(temp[0].header)
tempmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
colmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(tempmap[0].header)
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule

from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))

    freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3CN',
                                                             fmin=ch3cncube.spectral_axis.min(),
                                                             fmax=ch3cncube.spectral_axis.max())
    
    
    # TODO: determine this from data
    v_cen = 10*u.km/u.s
    v_disp = 1.5*u.km/u.s

    mod = lte_molecule.generate_model(ch3cncube.spectral_axis, v_cen, v_disp, temp, N_tot,
                              freqs, aij, deg, EU, partfunc)

    data_sp = ch3cncube[:, y, x]
    data_sp_K = data_sp.value * ch3cncube.jtok_factors()
    ax.plot(ch3cncube.spectral_axis.to(u.GHz), data_sp_K)
    ax.plot(ch3cncube.spectral_axis.to(u.GHz),
            mod, linestyle='-', color='k', linewidth=0.5)
########################################################
# Started Logging At: 2022-05-18 16:02:59
########################################################
########################################################
# # Started Logging At: 2022-05-18 16:02:59
########################################################
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
import numpy as np
# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
