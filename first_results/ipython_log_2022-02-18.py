########################################################
# Started Logging At: 2022-02-18 16:54:10
########################################################
########################################################
# # Started Logging At: 2022-02-18 16:54:10
########################################################
import time
start = time.time()

import pylab as pl
display_dpi = 150

# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
cubefns
#[Out]# ['/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_103_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_110_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_149_spw91_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_151_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_139_spw87_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_98_spw106_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_271_spw49_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_135_spw47_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_144_spw49_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_137_spw69_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_137_spw85_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_273_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_95_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_141_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_98_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_254_spw85_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_93_spw104_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_91_spw102_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_127_spw65_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_263_spw111_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_245_spw67_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_258_spw69_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_106_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_270_spw91_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_142_spw111_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_129_spw67_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_93_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_87_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_114_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_152_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_89_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_101_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_132_spw107_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_249_spw107_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_112_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_264_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_111_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_87_spw102_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_261_spw109_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw23_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_130_spw105_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_259_spw71_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_252_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_266_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_259_spw47_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw108_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_257_spw45_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_89_spw104_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_244_spw65_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_104_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_134_spw45_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_250_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw106_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_147_spw89_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_142_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_91_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_97_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_255_spw87_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_103_spw108_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_247_spw105_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_107_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_140_spw109_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_139_spw71_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_268_spw89_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image']
from astropy import coordinates, units as u
crd = coordinates.SkyCoord("17:46:10.59, -28:42:18.6", frame='fk5', unit=(u.h, u.deg))
crd = coordinates.SkyCoord("17:46:10.59 -28:42:18.6", frame='fk5', unit=(u.h, u.deg))
ch3cncube.wcs.celestial.world_to_pixel(crd)
for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    break
ch3cncube.wcs.celestial.world_to_pixel(crd)
#[Out]# (array(258.62486566), array(251.65342189))
get_ipython().run_line_magic('pinfo', 'ch3cncube.wcs.celestial.world_to_pixel')
get_ipython().run_line_magic('pinfo', 'ch3cncube.wcs.celestial.world_to_pixel_values')
for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
    pl.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
    pl.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
    break
for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        pl.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
len(cubes)
#[Out]# 8
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
from astropy import coordinates, units as u
from astropy import wcs
temp = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
col = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(temp[0].header)
from astropy import coordinates, units as u
from astropy import wcs
from astropy.io import fits
temp = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
col = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(temp[0].header)
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule

from spectral_cube import SpectralCube
from astropy import units as u
from lte_modeling_tools import get_molecular_parameters
from astropy import constants
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
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
    freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3CN', fmin=ch3cncube.spectral_axis.min(), fmax=ch3cn.spectral_axis.max())
    
    
    # TODO: determine this from data
    v_cen = 10*u.km/u.s
    v_disp = 1.5*u.km/u.s

    mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                              freqs, aij, deg, EU, partfunc)
    mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

    data_sp = ch3cncube[:, y, x]
    data_sp_K = data_sp.value * ch3cncube.jtok_factors()
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K.value, linestyle='-', color='k', linewidth=0.5)
tempmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
colmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(temp[0].header)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
    freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3CN', fmin=ch3cncube.spectral_axis.min(), fmax=ch3cn.spectral_axis.max())
    
    
    # TODO: determine this from data
    v_cen = 10*u.km/u.s
    v_disp = 1.5*u.km/u.s

    mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                              freqs, aij, deg, EU, partfunc)
    mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

    data_sp = ch3cncube[:, y, x]
    data_sp_K = data_sp.value * ch3cncube.jtok_factors()
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K.value, linestyle='-', color='k', linewidth=0.5)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
    freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3CN', fmin=ch3cncube.spectral_axis.min(), fmax=ch3cncube.spectral_axis.max())
    
    
    # TODO: determine this from data
    v_cen = 10*u.km/u.s
    v_disp = 1.5*u.km/u.s

    mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                              freqs, aij, deg, EU, partfunc)
    mod_sp = pyspeckit.Spectrum(xarr=sp_axis.to(u.GHz), data = mod, unit = u.K)

    data_sp = ch3cncube[:, y, x]
    data_sp_K = data_sp.value * ch3cncube.jtok_factors()
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K.value, linestyle='-', color='k', linewidth=0.5)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
    freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3CN',
                                                             fmin=ch3cncube.spectral_axis.min(),
                                                             fmax=ch3cncube.spectral_axis.max())
    
    
    # TODO: determine this from data
    v_cen = 10*u.km/u.s
    v_disp = 1.5*u.km/u.s

    mod = lte_molecule.generate_model(sp_axis, v_cen, v_disp, temp, N_tot,
                              freqs, aij, deg, EU, partfunc)

    data_sp = ch3cncube[:, y, x]
    data_sp_K = data_sp.value * ch3cncube.jtok_factors()
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K.value, linestyle='-', color='k', linewidth=0.5)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
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
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K.value, linestyle='-', color='k', linewidth=0.5)
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis, ch3cncube[:,y,x].value)
    
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
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K, linestyle='-', color='k', linewidth=0.5)
import pylab as pl
pl.rcParams['fig.facecolor'] = 'w'

# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'

# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis.to(u.GHz), ch3cncube[:,y,x].value)
    
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
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K, linestyle='-', color='k', linewidth=0.5)
    break
pl.figure(figsize=(12,12))
for ind, ch3cncube in enumerate(cubes):
    ax = pl.subplot(3,3,ind+1)
    
    x,y = map(int, ww.celestial.world_to_pixel(crd))
    T = temp = tempmap[0].data[y,x]
    N = N_tot = colmap[0].data[y,x]
        
    x,y = map(int, ch3cncube.wcs.celestial.world_to_pixel(crd))
        
    ax.plot(ch3cncube.spectral_axis.to(u.GHz), ch3cncube[:,y,x].value)
    
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
    pl.plot(ch3cncube.spectral_axis.to(u.GHz),
            data_sp_K, linestyle='-', color='k', linewidth=0.5)
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
    if ii == 2:
        break
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
    if ind == 2:
        break
mod
#[Out]# array([nan, nan, nan, ..., nan, nan, nan])
temp
#[Out]# 0.0
N
#[Out]# 0.0
crd = coordinates.SkyCoord("17:46:11.05 -28:42:14.6", frame='fk5', unit=(u.h, u.deg))
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
    if ind == 2:
        break
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
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
        print(ch3cntbl)
import warnings
warnings.filterwarnings('ignore')
N,T
#[Out]# (13.897920720743063, 495.52107936939143)
crd = coordinates.SkyCoord("17:46:10.80 -28:42:13.3", frame='fk5', unit=(u.h, u.deg))
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
N,T
#[Out]# (13.740148490539461, 153.10917256985616)
ch3cncube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#  n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#  n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#  n_s:   3840  type_s: FREQ      unit_s: Hz     range: 255308020805.109 Hz:257182502433.963 Hz
ch3cncube.jtok_factors()
#[Out]# array([11.708251  , 11.70825372, 11.70825   , ..., 11.70823204,
#[Out]#        11.70823613, 11.70828092])
cubefns
#[Out]# ['/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_103_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_110_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_149_spw91_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_151_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_139_spw87_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_98_spw106_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_271_spw49_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_135_spw47_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_144_spw49_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_137_spw69_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_137_spw85_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_273_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_95_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_141_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_98_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_254_spw85_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_93_spw104_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_91_spw102_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_127_spw65_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_263_spw111_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_245_spw67_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_258_spw69_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_106_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_270_spw91_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_142_spw111_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_129_spw67_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_93_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_146_spw51_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_87_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_114_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_152_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_89_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_101_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_132_spw107_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_249_spw107_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_112_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_264_spw29_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_111_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_87_spw102_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_261_spw109_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw23_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_130_spw105_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_259_spw71_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_252_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_266_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_259_spw47_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_99_spw108_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_257_spw45_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_89_spw104_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_244_spw65_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_104_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_134_spw45_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_250_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_102_spw106_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_147_spw89_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_142_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_91_spw25_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_97_spw27_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_255_spw87_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_103_spw108_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_247_spw105_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_107_spw31_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_140_spw109_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_139_spw71_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image',
#[Out]#  '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/source_ab_268_spw89_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image']
cubes = []

for fn in cubefns:

    ch3cncube = SpectralCube.read(fn)
    
    ch3cntbl = ch3cncube.find_lines(chemical_name='CH3CN', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    ch3cntbl = ch3cntbl[ch3cntbl['Quantum Number Code'] == 202]
    if len(ch3cntbl) > 0:
        cubes.append(ch3cncube)
        print(fn)
        print(ch3cntbl)
