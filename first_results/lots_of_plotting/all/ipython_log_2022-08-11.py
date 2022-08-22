########################################################
# Started Logging At: 2022-08-11 13:29:30
########################################################
########################################################
# # Started Logging At: 2022-08-11 13:29:30
########################################################
# For parallelization, which helps to make convolution faster
import dask
dask.config.set(scheduler = 'threads', num_workers = 8)
from dask.diagnostics import ProgressBar
ProgressBar().register()
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u

results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# Get names of spectral windows
# I think these are in vaguely the right order?
freq_spws = [
'87_spw25',
# '87_spw102', # duplicate of above
'89_spw27',
# '89_spw104', # duplicate of above

'91_spw25', 
# '91_spw102', # duplicate of above
'93_spw27',
# '93_spw104', # duplicate of above
'95_spw25', 
'97_spw27', 
'98_spw29',
# '98_spw106', # duplicate of above
# '99_spw25', # duplicate of above
'99_spw31',  
# '99_spw108', # duplicate of above

# '101_spw27', # duplicate of above
'102_spw23',
'102_spw29', 
# '102_spw106', # duplicate of above
'104_spw25', # mislabeled, should be before 103_spw31
'103_spw31', 
# '103_spw108', # duplicate of above
'106_spw29',
'107_spw31',

'110_spw29',
'111_spw31',
'112_spw27',
'114_spw29',

'127_spw65', 
'129_spw67',

'130_spw105',
'132_spw107',
'134_spw45', 
'135_spw47',
'137_spw85',
'137_spw69',
# '139_spw87', # duplicate of above
'139_spw71',
# '141_spw25', # duplicate of above

'140_spw109',
'142_spw111',
'144_spw49',
'146_spw51',
'147_spw89',
'149_spw91',

'142_spw27', # mislabeled
# '151_spw29', # duplicate of above
'152_spw31',

'244_spw65',
'245_spw67',
'247_spw105',
'249_spw107',

'250_spw25',
'252_spw27',
'254_spw85',
'255_spw87',
'257_spw45',
'259_spw47',
# '258_spw69', # significant overlap with above?
'259_spw71',

'261_spw109',
'263_spw111',
'264_spw29',
'266_spw31',
'268_spw89',

'270_spw91',
'271_spw49',
'273_spw51'
]
# Only look at one cube from here on out
ind = 0
from spectral_cube import SpectralCube

cube = SpectralCube.read(results+'source_ab_'+freq_spws[ind]+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image', use_dask = True)
# Extract spectrum for a coordinate in central core region
from astropy import coordinates
from astropy import wcs

crd = coordinates.SkyCoord("17:46:10.63 -28:42:17.8", frame='fk5', unit=(u.h, u.deg))

x, y = map(int, smoothed_cube.wcs.celestial.world_to_pixel(crd))
########################################################
# Started Logging At: 2022-08-11 13:51:07
########################################################
########################################################
# # Started Logging At: 2022-08-11 13:51:07
########################################################
# For parallelization, which helps to make convolution faster
import dask
dask.config.set(scheduler = 'threads', num_workers = 8)
from dask.diagnostics import ProgressBar
ProgressBar().register()
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u

results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# Get names of spectral windows
# I think these are in vaguely the right order?
freq_spws = [
'87_spw25',
# '87_spw102', # duplicate of above
'89_spw27',
# '89_spw104', # duplicate of above

'91_spw25', 
# '91_spw102', # duplicate of above
'93_spw27',
# '93_spw104', # duplicate of above
'95_spw25', 
'97_spw27', 
'98_spw29',
# '98_spw106', # duplicate of above
# '99_spw25', # duplicate of above
'99_spw31',  
# '99_spw108', # duplicate of above

# '101_spw27', # duplicate of above
'102_spw23',
'102_spw29', 
# '102_spw106', # duplicate of above
'104_spw25', # mislabeled, should be before 103_spw31
'103_spw31', 
# '103_spw108', # duplicate of above
'106_spw29',
'107_spw31',

'110_spw29',
'111_spw31',
'112_spw27',
'114_spw29',

'127_spw65', 
'129_spw67',

'130_spw105',
'132_spw107',
'134_spw45', 
'135_spw47',
'137_spw85',
'137_spw69',
# '139_spw87', # duplicate of above
'139_spw71',
# '141_spw25', # duplicate of above

'140_spw109',
'142_spw111',
'144_spw49',
'146_spw51',
'147_spw89',
'149_spw91',

'142_spw27', # mislabeled
# '151_spw29', # duplicate of above
'152_spw31',

'244_spw65',
'245_spw67',
'247_spw105',
'249_spw107',

'250_spw25',
'252_spw27',
'254_spw85',
'255_spw87',
'257_spw45',
'259_spw47',
# '258_spw69', # significant overlap with above?
'259_spw71',

'261_spw109',
'263_spw111',
'264_spw29',
'266_spw31',
'268_spw89',

'270_spw91',
'271_spw49',
'273_spw51'
]
# Only look at one cube from here on out
ind = 0
from spectral_cube import SpectralCube

cube = SpectralCube.read(results+'source_ab_'+freq_spws[ind]+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image', use_dask = True)
cube_common_beam = cube.beams.common_beam(max_iter = 20, max_epsilon = 0.01)
smoothed_cube = cube.convolve_to(cube_common_beam) # Convert from VaryingResolution to regular

# import radio_beam

# beams = []
# for cube in smoothed_cubes:
#     beams.append(cube.beam)

# common_beam = radio_beam.commonbeam.common_manybeams_mve(radio_beam.Beams(beams = beams)) # Get common beam between all cubes
# common_beam
# resampled_cubes = []
# for cube in smoothed_cubes:
#     resampled_cube = cube.convolve_to(common_beam) # This cell takes almost no time; the cubes do not get convolved in this step
#     resampled_cubes.append(resampled_cube)
# Extract spectrum for a coordinate in central core region
from astropy import coordinates
from astropy import wcs

crd = coordinates.SkyCoord("17:46:10.63 -28:42:17.8", frame='fk5', unit=(u.h, u.deg))

x, y = map(int, smoothed_cube.wcs.celestial.world_to_pixel(crd))
# # Extract spectrum for brightest pixel in central core region
# import regions
# reg = regions.Regions.read('/blue/adamginsburg/abulatek/brick/first_results/region_temperatures/centralcoreellipse.reg')

# smoothed_cube_reg = smoothed_cube.subcube_from_regions(reg)
# m0 = smoothed_cube.moment0()
# Take the mean spectrum of the outflow cavity region
import regions

reg = regions.Regions.read('outflow_cavity_region.reg')
smoothed_scube = smoothed_cube.subcube_from_regions(reg)
mean_spectrum = smoothed_scube.mean(axis = (1, 2)) #how='slice', progressbar=True
# For parallelization, which helps to make convolution faster
import dask
dask.config.set(scheduler = 'threads', num_workers = 8)
from dask.diagnostics import ProgressBar
ProgressBar().register()
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u

results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
# Get names of spectral windows
# I think these are in vaguely the right order?
freq_spws = [
'87_spw25',
# '87_spw102', # duplicate of above
'89_spw27',
# '89_spw104', # duplicate of above

'91_spw25', 
# '91_spw102', # duplicate of above
'93_spw27',
# '93_spw104', # duplicate of above
'95_spw25', 
'97_spw27', 
'98_spw29',
# '98_spw106', # duplicate of above
# '99_spw25', # duplicate of above
'99_spw31',  
# '99_spw108', # duplicate of above

# '101_spw27', # duplicate of above
'102_spw23',
'102_spw29', 
# '102_spw106', # duplicate of above
'104_spw25', # mislabeled, should be before 103_spw31
'103_spw31', 
# '103_spw108', # duplicate of above
'106_spw29',
'107_spw31',

'110_spw29',
'111_spw31',
'112_spw27',
'114_spw29',

'127_spw65', 
'129_spw67',

'130_spw105',
'132_spw107',
'134_spw45', 
'135_spw47',
'137_spw85',
'137_spw69',
# '139_spw87', # duplicate of above
'139_spw71',
# '141_spw25', # duplicate of above

'140_spw109',
'142_spw111',
'144_spw49',
'146_spw51',
'147_spw89',
'149_spw91',

'142_spw27', # mislabeled
# '151_spw29', # duplicate of above
'152_spw31',

'244_spw65',
'245_spw67',
'247_spw105',
'249_spw107',

'250_spw25',
'252_spw27',
'254_spw85',
'255_spw87',
'257_spw45',
'259_spw47',
# '258_spw69', # significant overlap with above?
'259_spw71',

'261_spw109',
'263_spw111',
'264_spw29',
'266_spw31',
'268_spw89',

'270_spw91',
'271_spw49',
'273_spw51'
]
# Only look at one cube from here on out
ind = freq_spws.index('129_spw67')
from spectral_cube import SpectralCube

cube = SpectralCube.read(results+'source_ab_'+freq_spws[ind]+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image', use_dask = True)
cube_common_beam = cube.beams.common_beam(max_iter = 20, max_epsilon = 0.01)
smoothed_cube = cube.convolve_to(cube_common_beam) # Convert from VaryingResolution to regular

# import radio_beam

# beams = []
# for cube in smoothed_cubes:
#     beams.append(cube.beam)

# common_beam = radio_beam.commonbeam.common_manybeams_mve(radio_beam.Beams(beams = beams)) # Get common beam between all cubes
# common_beam
# resampled_cubes = []
# for cube in smoothed_cubes:
#     resampled_cube = cube.convolve_to(common_beam) # This cell takes almost no time; the cubes do not get convolved in this step
#     resampled_cubes.append(resampled_cube)
# Extract spectrum for a coordinate in central core region
from astropy import coordinates
from astropy import wcs

crd = coordinates.SkyCoord("17:46:10.63 -28:42:17.8", frame='fk5', unit=(u.h, u.deg))

x, y = map(int, smoothed_cube.wcs.celestial.world_to_pixel(crd))
# # Extract spectrum for brightest pixel in central core region
# import regions
# reg = regions.Regions.read('/blue/adamginsburg/abulatek/brick/first_results/region_temperatures/centralcoreellipse.reg')

# smoothed_cube_reg = smoothed_cube.subcube_from_regions(reg)
# m0 = smoothed_cube.moment0()
# Take the mean spectrum of the outflow cavity region
import regions

reg = regions.Regions.read('outflow_cavity_region.reg')
smoothed_scube = smoothed_cube.subcube_from_regions(reg)
#notused mean_spectrum = smoothed_scube.mean(axis = (1, 2)) #how='slice', progressbar=True
import matplotlib.pyplot as plt

SM_SIZE = 12
MD_SIZE = 18
LG_SIZE = 22

plt.rc('font', size = MD_SIZE)          # controls default text sizes
plt.rc('axes', titlesize = LG_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize = MD_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize = MD_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize = MD_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize = SM_SIZE)    # legend fontsize
plt.rc('figure', titlesize = LG_SIZE)   # fontsize of the figure title
ax = plt.subplot(projection = m0.wcs)
plt.plot(x, y, marker = 'x', color = 'r')
im = ax.imshow(m0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label(f'Integrated Intensity [{m0.unit}]')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
m0 = smoothed_cube.moment0()
#m0 = smoothed_cube.moment0()
m0 = smoothed_cube[0] # HACK
ax = plt.subplot(projection = m0.wcs)
plt.plot(x, y, marker = 'x', color = 'r')
im = ax.imshow(m0.value, origin='lower', cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label(f'Integrated Intensity [{m0.unit}]')
ax.set_ylabel('Declination')
ax.set_xlabel('Right Ascension')
spectrum = smoothed_cube[:, y, x].to(u.K)
import numpy as np

fig = plt.figure(figsize = (30, 8))

plt.plot(smoothed_cube.spectral_axis.to(u.GHz), np.array(spectrum), linestyle = '-', color = 'k', linewidth = 1, drawstyle = 'steps-mid')
plt.xlabel(f"Frequency [{smoothed_cube.spectral_axis.to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
plt.xlim(smoothed_cube.spectral_extrema.to(u.GHz).value)
plt.title(f"{freq_spws[ind]} data spectrum")
plt.show()
# Get path of XCLASS directory
import os
import sys
XCLASSRootDir = str(os.environ.get('XCLASSRootDir', '')).strip()
XCLASSRootDir = os.path.normpath(XCLASSRootDir) + "/"

# Extend sys.path variable
NewPath = XCLASSRootDir + "build_tasks/"
if (not NewPath in sys.path):
    sys.path.append(NewPath)

import task_myXCLASS
# source_size =  0.33 # total guess, but based on Dan's data? also, not used
tkin = 151 # 7.9842254347e+01
Ntot = 5.0e+13 # 10**N_tot # 1.5061155735e+14 # 1.0061155735e+16
# vwidth_fwhm = 4.21
# vwidth = vwidth_fwhm / np.sqrt(8*np.log(2))
# vcen = 38.0
# tbg = 0
# Define path and name of molfit file
# DEFAULT MolfitsFileName = LocalPath + "files/my_molecules.molfit"
LocalPath = os.getcwd() + "/"

MolfitsFileName = LocalPath + 'my_molecules-tout.molfit'

# def set_xclass_parameters(source_size=source_size, tkin=tkin, Ntot=Ntot, vwidth_fwhm=vwidth_fwhm, vcen=vcen, filename=MolfitsFileName):
#     with open(filename, "w") as fh:
#         fh.write(f"""
# %============================================================================================================================
# %
# % define molecules and their components:
# %
# %============================================================================================================================
# %    source size:                 T_kin:               N_tot:            V_width:                V_off:   CFFlag:    keyword:
# %       [arcsec]:                   [K]:              [cm-2]:            [km /s]:              [km /s]:       []:         []:
# CH3CN;v=0;           1
#  {source_size}       {tkin}     {Ntot}     {vwidth_fwhm}    {vcen}         c
# """)
# define freq. step (in MHz)
FreqStep = 0.1

# define beam minor axis length (in arsec)
BMIN = smoothed_cube.beam.minor.to(u.arcsec).value # None

# define beam major axis length (in arsec)
BMAJ = smoothed_cube.beam.major.to(u.arcsec).value # None

# define beam position angle (in degree)
BPA = smoothed_cube.beam.pa.to(u.deg).value # None

# depending on parameter "Inter_Flag" define beam size (in arcsec)
# (Inter_Flag = True) or size of telescope (in m) (Inter_Flag = False)
TelescopeSize = np.abs(BMIN**2 - BMAJ**2)**0.5 # arcsec # 1000.0 # meters

# interferometric data?
Inter_Flag = True

# define red shift
Redshift = None

# BACKGROUND: describe continuum with tBack and tslope only
t_back_flag = True

# BACKGROUND: define background temperature (in K)
tBack = 0.0

# BACKGROUND: define temperature slope (dimensionless)
tslope = 0.0

# BACKGROUND: define path and name of ASCII file describing continuum as function
#             of frequency
BackgroundFileName = ""

# DUST: define hydrogen column density (in cm^(-2))
N_H = 1.e22

# DUST: define spectral index for dust (dimensionless)
beta_dust = 0.0

# DUST: define kappa at 1.3 mm (cm^(2) g^(-1))
kappa_1300 = 0.0

# DUST: define path and name of ASCII file describing dust opacity as
#       function of frequency
DustFileName = ""

# FREE-FREE: define electronic temperature (in K)
Te_ff = None

# FREE-FREE: define emission measure (in pc cm^(-6))
EM_ff = None

# SYNCHROTRON: define kappa of energy spectrum of electrons (electrons m^(−3) GeV^(-1))
kappa_sync = None

# SYNCHROTRON: define magnetic field (in Gauss)
B_sync = None

# SYNCHROTRON: energy spectral index (dimensionless)
p_sync = None

# SYNCHROTRON: thickness of slab (in AU)
l_sync = None

# PHEN-CONT: define phenomenological function which is used to describe
#            the continuum
ContPhenFuncID = None

# PHEN-CONT: define first parameter for phenomenological function
ContPhenFuncParam1 = None

# PHEN-CONT: define second parameter for phenomenological function
ContPhenFuncParam2 = None

# PHEN-CONT: define third parameter for phenomenological function
ContPhenFuncParam3 = None

# PHEN-CONT: define fourth parameter for phenomenological function
ContPhenFuncParam4 = None

# PHEN-CONT: define fifth parameter for phenomenological function
ContPhenFuncParam5 = None

# use iso ratio file?
iso_flag = True

# define path and name of iso ratio file
#DEFAULT IsoTableFileName = LocalPath + "files/my_isonames.txt"
IsoTableFileName = LocalPath + "my_isonames.txt"

# define path and name of file describing Non-LTE parameters
CollisionFileName = ""

# define number of pixels in x-direction (used for sub-beam description)
NumModelPixelXX = 100

# define number of pixels in y-direction (used for sub-beam description)
NumModelPixelYY = 100

# take local-overlap into account or not
LocalOverlapFlag = False

# disable sub-beam description
NoSubBeamFlag = True

# define path and name of database file
dbFilename = ""

# define rest freq. (in MHz)
RestFreq = 0.0

# define v_lsr (in km/s)
vLSR = 0.0
import io
from contextlib import redirect_stdout

def myxclass_call(FreqMin=1e3, FreqMax=1e4, verbose=False):
    ## call myXCLASS function
    with io.StringIO() as buf, redirect_stdout(buf):
        modeldata, log, TransEnergies, IntOpt, JobDir = task_myXCLASS.myXCLASS(
                                                    FreqMin, FreqMax, FreqStep,
                                                    TelescopeSize, BMIN, BMAJ,
                                                    BPA, Inter_Flag, Redshift,
                                                    t_back_flag, tBack, tslope,
                                                    BackgroundFileName,
                                                    N_H, beta_dust, kappa_1300,
                                                    DustFileName, Te_ff, EM_ff,
                                                    kappa_sync, B_sync, p_sync,
                                                    l_sync, ContPhenFuncID,
                                                    ContPhenFuncParam1,
                                                    ContPhenFuncParam2,
                                                    ContPhenFuncParam3,
                                                    ContPhenFuncParam4,
                                                    ContPhenFuncParam5,
                                                    MolfitsFileName, iso_flag,
                                                    IsoTableFileName,
                                                    CollisionFileName,
                                                    NumModelPixelXX,
                                                    NumModelPixelYY,
                                                    LocalOverlapFlag,
                                                    NoSubBeamFlag,
                                                    dbFilename,
                                                    RestFreq, vLSR)
        output = buf.getvalue()
    if verbose:
        print(output)
        
    return modeldata, log, TransEnergies, IntOpt, JobDir
# define min. freq. (in MHz)
FreqMin = np.min(spectrum.spectral_axis.to(u.MHz)).value
# define max. freq. (in MHz)
FreqMax = np.max(spectrum.spectral_axis.to(u.MHz)).value
# # set_xclass_parameters() # The molfit file is already created, I don't want to mess with it; to run this, need to uncomment function definition a few cells ago
# modeldata, log, TransEnergies, IntOpt, JobDir = myxclass_call(FreqMin=FreqMin, FreqMax=FreqMax)
# fig = plt.figure(figsize = (30, 8))

# plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], linestyle = '-', color = 'r', linewidth = 1, drawstyle = 'steps-mid')
# plt.xlabel(f"Frequency [{(modeldata[:,0]*u.MHz).to(u.GHz).unit}]")
# plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
# plt.xlim(smoothed_cube.spectral_extrema.to(u.GHz).value)
# plt.title(f"{freq_spws[ind]} model spectrum, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92
# plt.show()
# Preliminary continuum subtraction
spectrum_contsub = np.array(spectrum) - np.median(np.array(spectrum))
# fig = plt.figure(figsize = (30, 8))

# plt.plot(smoothed_cube.spectral_axis.to(u.GHz), spectrum_contsub, linestyle = '-', color = 'k', linewidth = 0.75, drawstyle = 'steps-mid', 
#          label = "Data")
# plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], linestyle = '--', color = 'r', linewidth = 1, drawstyle = 'steps-mid', 
#          label = "Model")
# plt.legend(loc='best')
# plt.xlabel(f"Frequency [{smoothed_cube.spectral_axis.to(u.GHz).unit}]")
# plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
# plt.xlim(smoothed_cube.spectral_extrema.to(u.GHz).value)
# plt.title(f"{freq_spws[ind]} data and model comparison, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92
# plt.show()
import glob
MolfitsFileNames = glob.glob(LocalPath + 'individuals/*.molfit')
modeldata_all = []
transenergies_all = []
for MolfitsFileName in MolfitsFileNames:
    print(MolfitsFileName)
    modeldata, log, TransEnergies, IntOpt, JobDir = myxclass_call(FreqMin=FreqMin, FreqMax=FreqMax)
    modeldata_all.append(modeldata)
    transenergies_all.append(TransEnergies)
# User-maintained color list
colors = {'13CO;v=0;'   : ['tab:blue',   '--'],
          'CCCS;v=0;'   : ['tab:orange', '--'],
          'CCH;v=0;'    : ['tab:green',  '--'],
          'CCS;v=0;'    : ['tab:red',    '--'],
          'CH3CN;v=0;'  : ['tab:purple', '--'],
          'CH3OH;v=0;'  : ['tab:pink',   '--'],
          'CO;v=0;'     : ['tab:olive',  '--'],
          'H2CS;v=0;'   : ['tab:cyan',   '--'],
          'HC3N;v=0;'   : ['tab:blue',   '-.'],
          'HCN;v=0;'    : ['tab:orange', '-.'],
          'HCO+;v=0;'   : ['tab:green',  '-.'],
          'HDCO;v=0;'   : ['tab:red',    '-.'],
          'HNCO;v=0;'   : ['tab:purple', '-.'],
          'HNC;v=0;'    : ['tab:pink',   '-.'],
          'HOCO+;v=0;'  : ['tab:olive',  '-.'],
          'N2H+;v=0;'   : ['tab:cyan',   '-.'],
          'NH2D;v=0;'   : ['tab:blue',   ':'],
          'NH3;v=0;'    : ['tab:orange', ':'],
          'OCS;v=0;'    : ['tab:green',  ':'],
          'C17O;v=0;'   : ['tab:red',    ':'],
          'H13CO+;v=0;' : ['tab:purple', ':'],
          '13CO;v=0;'   : ['tab:pink',   ':'],
          'HC15N;v=0;'  : ['tab:olive',  ':'],
          'HN13C;v=0;'  : ['tab:cyan',   ':'],
          'C18O;v=0;'   : ['tab:blue',   (0,(3,1,1,1,1,1))],
          #'H2COH+;v=0;' : ['tab:orange', (0,(3,1,1,1,1,1))],
          'NH2CHO;v=0;' : ['tab:green',  (0,(3,1,1,1,1,1))],
          'H13CN;v=0;'  : ['tab:red',    (0,(3,1,1,1,1,1))],
         }
# to add: C3H2, CH3OCHO, CH3OCH3, H2CO, SO, SO2, 13CS, CS, CN
# column densities: {'CO': 1e18, '13CO': 1e17, 'H2CO': 1e16, 'CH3OH': 1e15, 'CH3OCHO': 1e14, 'CH3OCH3': }
# fig = plt.figure(figsize = (30, 8))

# for index, (modeldata, transenergies) in enumerate(zip(modeldata_all, transenergies_all)):
#     if len(modeldata) != 0:
#         colorname = transenergies[1][-1]
#         molname = MolfitsFileNames[index].replace('/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/', '').replace('.molfit', '')
#         plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], colors[colorname], linewidth = 1, drawstyle = 'steps-mid', 
#                  label = molname)
#         plt.xlabel(f"Frequency [{(modeldata[:,0]*u.MHz).to(u.GHz).unit}]")
#     else:
#         print(f"***** Warning! {molname} was not modeled! *****")
# plt.legend(loc = 'best')
# plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
# plt.xlim(smoothed_cube.spectral_extrema.to(u.GHz).value)
# plt.title(f"{freq_spws[ind]} model spectrum, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92 # NOTE: TKIN AND NTOT MAY CHANGE WHEN WE FIT !
# plt.show()
from astropy.stats import mad_std

# Get MAD-estimated RMS of data spectrum
mad_rms = mad_std(spectrum_contsub)
# # Set up color and linestyle cycler
# from cycler import cycler
# plt.rc('axes', prop_cycle=(cycler('color', ['xkcd:azure', 'xkcd:orange', 'xkcd:teal', 'xkcd:magenta']) * cycler('linestyle', ['--', ':', '-.', (0,(3,1,1,1,1,1))])))
fig = plt.figure(figsize = (30, 8))

plt.plot(smoothed_cube.spectral_axis.to(u.GHz), spectrum_contsub, linestyle = '-', color = 'grey', linewidth = 1, drawstyle = 'steps-mid', 
        label = "Data")
plt.axhspan(ymin = -mad_rms, ymax = mad_rms, alpha = 0.25, color = 'k', label = 'RMS level in data')
for index, (modeldata, transenergies) in enumerate(zip(modeldata_all, transenergies_all)):
    molname = MolfitsFileNames[index].replace('/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/', '').replace('.molfit', '')
    if len(modeldata) != 0:
        colorname = transenergies[1][-1]
        plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], color = colors[colorname][0], linestyle = colors[colorname][1], linewidth = 1, 
                 drawstyle = 'steps-mid', label = molname)
    else:
        print(f"***** Warning! {molname} was not modeled! *****")
plt.legend(loc = 'best')
plt.xlabel(f"Frequency [{smoothed_cube.spectral_axis.to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
plt.xlim(smoothed_cube.spectral_extrema.to(u.GHz).value)
# plt.xlim(137.25, 137.50)
plt.title(f"{freq_spws[ind]} data and model comparison, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92 # NOTE: TKIN AND NTOT MAY CHANGE WHEN WE FIT !
plt.savefig(f"spectra_ID/{freq_spws[ind]}.pdf", facecolor = 'w', edgecolor = 'w', bbox_inches = 'tight') # dpi = 300
plt.show()
FreqMin, FreqMax
#[Out]# (126732.91597266681, 128607.48200316788)
import io
from contextlib import redirect_stdout

def myxclass_call(FreqMin=1e3, FreqMax=1e4, verbose=False):
    ## call myXCLASS function
    print(MolfitsFileName)
    with io.StringIO() as buf, redirect_stdout(buf):
        modeldata, log, TransEnergies, IntOpt, JobDir = task_myXCLASS.myXCLASS(
                                                    FreqMin, FreqMax, FreqStep,
                                                    TelescopeSize, BMIN, BMAJ,
                                                    BPA, Inter_Flag, Redshift,
                                                    t_back_flag, tBack, tslope,
                                                    BackgroundFileName,
                                                    N_H, beta_dust, kappa_1300,
                                                    DustFileName, Te_ff, EM_ff,
                                                    kappa_sync, B_sync, p_sync,
                                                    l_sync, ContPhenFuncID,
                                                    ContPhenFuncParam1,
                                                    ContPhenFuncParam2,
                                                    ContPhenFuncParam3,
                                                    ContPhenFuncParam4,
                                                    ContPhenFuncParam5,
                                                    MolfitsFileName, iso_flag,
                                                    IsoTableFileName,
                                                    CollisionFileName,
                                                    NumModelPixelXX,
                                                    NumModelPixelYY,
                                                    LocalOverlapFlag,
                                                    NoSubBeamFlag,
                                                    dbFilename,
                                                    RestFreq, vLSR)
        output = buf.getvalue()
    if verbose:
        print(output)
        
    return modeldata, log, TransEnergies, IntOpt, JobDir
MolfitsFileNames
#[Out]# ['/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HCNv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HDCOv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HCO+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/13COv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/NH3v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HNCv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/NH2Dv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HNCOv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/N2H+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCHv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HOCO+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CH3CNv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HC3Nv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/COv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CH3OHv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/OCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/H2CSv=0.molfit']
MolfitsFileNames[12]
#[Out]# '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CH3CNv=0.molfit'
MolfitsFileNames
#[Out]# ['/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HCNv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HDCOv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HCO+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/13COv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/NH3v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HNCv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/NH2Dv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HNCOv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/N2H+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCHv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HOCO+v=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CH3CNv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HC3Nv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/COv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/CH3OHv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/OCSv=0.molfit',
#[Out]#  '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/H2CSv=0.molfit']
MolfitsFileNames[13]
#[Out]# '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HC3Nv=0.molfit'
MolfitsFileName = '/blue/adamginsburg/abulatek/brick/first_results/lots_of_plotting/all/individuals/HC3Nv=0.molfit'
modeldata, log, TransEnergies, IntOpt, JobDir = myxclass_call(FreqMin=FreqMin, FreqMax=FreqMax)
# Take the mean spectrum of the outflow cavity region
import regions

reg = regions.Regions.read('outflow_cavity_region.reg')
smoothed_scube = smoothed_cube.subcube_from_regions(reg)
mean_spectrum = smoothed_scube.mean(axis = (1, 2)) #how='slice', progressbar=True
# THIS IS ADAM'S HACKERY
import pyspeckit
sp = pyspeckit.Spectrum(xarr=mean_spectrum.spectral_axis, data=mean_spectrum)
sp.plotter()
sp.plotter(xmin=127.26*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[1, 127.36766, 0.001])
sp.plotter(xmin=127.28*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[1, 127.36766, 0.001])
sp.plotter(xmin=127.29*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[1, 127.36766, 0.001])
sp.plotter(xmin=127.3*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[1, 127.36766, 0.001])
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[1, 127.36766, 0.001])
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
sp.xarr.convert_to(u.GHz)
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
sp.xarr.convert(u.GHz)
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
dir(sp.xarr)
#[Out]# ['T',
#[Out]#  '__abs__',
#[Out]#  '__add__',
#[Out]#  '__and__',
#[Out]#  '__array__',
#[Out]#  '__array_finalize__',
#[Out]#  '__array_function__',
#[Out]#  '__array_interface__',
#[Out]#  '__array_prepare__',
#[Out]#  '__array_priority__',
#[Out]#  '__array_struct__',
#[Out]#  '__array_ufunc__',
#[Out]#  '__array_wrap__',
#[Out]#  '__bool__',
#[Out]#  '__class__',
#[Out]#  '__class_getitem__',
#[Out]#  '__complex__',
#[Out]#  '__contains__',
#[Out]#  '__copy__',
#[Out]#  '__deepcopy__',
#[Out]#  '__delattr__',
#[Out]#  '__delitem__',
#[Out]#  '__dict__',
#[Out]#  '__dir__',
#[Out]#  '__divmod__',
#[Out]#  '__dlpack__',
#[Out]#  '__dlpack_device__',
#[Out]#  '__doc__',
#[Out]#  '__eq__',
#[Out]#  '__float__',
#[Out]#  '__floordiv__',
#[Out]#  '__format__',
#[Out]#  '__ge__',
#[Out]#  '__getattr__',
#[Out]#  '__getattribute__',
#[Out]#  '__getitem__',
#[Out]#  '__gt__',
#[Out]#  '__hash__',
#[Out]#  '__iadd__',
#[Out]#  '__iand__',
#[Out]#  '__ifloordiv__',
#[Out]#  '__ilshift__',
#[Out]#  '__imatmul__',
#[Out]#  '__imod__',
#[Out]#  '__imul__',
#[Out]#  '__index__',
#[Out]#  '__init__',
#[Out]#  '__init_subclass__',
#[Out]#  '__int__',
#[Out]#  '__invert__',
#[Out]#  '__ior__',
#[Out]#  '__ipow__',
#[Out]#  '__irshift__',
#[Out]#  '__isub__',
#[Out]#  '__iter__',
#[Out]#  '__itruediv__',
#[Out]#  '__ixor__',
#[Out]#  '__le__',
#[Out]#  '__len__',
#[Out]#  '__lshift__',
#[Out]#  '__lt__',
#[Out]#  '__matmul__',
#[Out]#  '__mod__',
#[Out]#  '__module__',
#[Out]#  '__mul__',
#[Out]#  '__ne__',
#[Out]#  '__neg__',
#[Out]#  '__new__',
#[Out]#  '__or__',
#[Out]#  '__pos__',
#[Out]#  '__pow__',
#[Out]#  '__quantity_subclass__',
#[Out]#  '__radd__',
#[Out]#  '__rand__',
#[Out]#  '__rdivmod__',
#[Out]#  '__reduce__',
#[Out]#  '__reduce_ex__',
#[Out]#  '__repr__',
#[Out]#  '__rfloordiv__',
#[Out]#  '__rlshift__',
#[Out]#  '__rmatmul__',
#[Out]#  '__rmod__',
#[Out]#  '__rmul__',
#[Out]#  '__ror__',
#[Out]#  '__rpow__',
#[Out]#  '__rrshift__',
#[Out]#  '__rshift__',
#[Out]#  '__rsub__',
#[Out]#  '__rtruediv__',
#[Out]#  '__rxor__',
#[Out]#  '__setattr__',
#[Out]#  '__setitem__',
#[Out]#  '__setstate__',
#[Out]#  '__sizeof__',
#[Out]#  '__str__',
#[Out]#  '__sub__',
#[Out]#  '__subclasshook__',
#[Out]#  '__truediv__',
#[Out]#  '__xor__',
#[Out]#  '_decompose',
#[Out]#  '_default_unit',
#[Out]#  '_dxarr',
#[Out]#  '_equivalencies',
#[Out]#  '_include_easy_conversion_members',
#[Out]#  '_make_header',
#[Out]#  '_new_view',
#[Out]#  '_not_implemented_or_raise',
#[Out]#  '_recursively_apply',
#[Out]#  '_repr_latex_',
#[Out]#  '_result_as_quantity',
#[Out]#  '_set_unit',
#[Out]#  '_to_own_unit',
#[Out]#  '_to_value',
#[Out]#  '_unit',
#[Out]#  '_unitstr',
#[Out]#  '_wrap_function',
#[Out]#  'all',
#[Out]#  'any',
#[Out]#  'argmax',
#[Out]#  'argmin',
#[Out]#  'argpartition',
#[Out]#  'argsort',
#[Out]#  'as_unit',
#[Out]#  'astype',
#[Out]#  'base',
#[Out]#  'byteswap',
#[Out]#  'cdelt',
#[Out]#  'center_frequency',
#[Out]#  'center_frequency_unit',
#[Out]#  'cgs',
#[Out]#  'choose',
#[Out]#  'clip',
#[Out]#  'compress',
#[Out]#  'conj',
#[Out]#  'conjugate',
#[Out]#  'convert_to_unit',
#[Out]#  'coord_to_x',
#[Out]#  'copy',
#[Out]#  'ctypes',
#[Out]#  'cumprod',
#[Out]#  'cumsum',
#[Out]#  'data',
#[Out]#  'decompose',
#[Out]#  'diagonal',
#[Out]#  'diff',
#[Out]#  'dot',
#[Out]#  'dtype',
#[Out]#  'dump',
#[Out]#  'dumps',
#[Out]#  'dxarr',
#[Out]#  'ediff1d',
#[Out]#  'equivalencies',
#[Out]#  'fill',
#[Out]#  'find_equivalencies',
#[Out]#  'flags',
#[Out]#  'flat',
#[Out]#  'flatten',
#[Out]#  'frame',
#[Out]#  'getfield',
#[Out]#  'imag',
#[Out]#  'in_range',
#[Out]#  'info',
#[Out]#  'insert',
#[Out]#  'isscalar',
#[Out]#  'item',
#[Out]#  'itemset',
#[Out]#  'itemsize',
#[Out]#  'make_dxarr',
#[Out]#  'max',
#[Out]#  'mean',
#[Out]#  'min',
#[Out]#  'nansum',
#[Out]#  'nbytes',
#[Out]#  'ndim',
#[Out]#  'newbyteorder',
#[Out]#  'nonzero',
#[Out]#  'partition',
#[Out]#  'prod',
#[Out]#  'ptp',
#[Out]#  'put',
#[Out]#  'ravel',
#[Out]#  'real',
#[Out]#  'redshift',
#[Out]#  'refX',
#[Out]#  'refX_unit',
#[Out]#  'refX_units',
#[Out]#  'repeat',
#[Out]#  'reshape',
#[Out]#  'resize',
#[Out]#  'round',
#[Out]#  'searchsorted',
#[Out]#  'set_unit',
#[Out]#  'setfield',
#[Out]#  'setflags',
#[Out]#  'shape',
#[Out]#  'si',
#[Out]#  'size',
#[Out]#  'sort',
#[Out]#  'squeeze',
#[Out]#  'std',
#[Out]#  'strides',
#[Out]#  'sum',
#[Out]#  'swapaxes',
#[Out]#  'take',
#[Out]#  'to',
#[Out]#  'to_string',
#[Out]#  'to_value',
#[Out]#  'tobytes',
#[Out]#  'tofile',
#[Out]#  'tolist',
#[Out]#  'tostring',
#[Out]#  'trace',
#[Out]#  'transpose',
#[Out]#  'umax',
#[Out]#  'umin',
#[Out]#  'unit',
#[Out]#  'units',
#[Out]#  'validate_unit',
#[Out]#  'value',
#[Out]#  'var',
#[Out]#  'velocity_convention',
#[Out]#  'view',
#[Out]#  'wcshead',
#[Out]#  'x_to_coord',
#[Out]#  'x_to_pix',
#[Out]#  'xtype']
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.32*u.GHz, xmax=127.4*u.GHz)
sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
sp.specfit.parinfo
#[Out]# [Param #0   AMPLITUDE0 =    0.0272988 +/-     0.000653408 ,
#[Out]#  Param #1       SHIFT0 =      127.351 +/-     2.78306e-05 ,
#[Out]#  Param #2       WIDTH0 =   0.00100696 +/-      2.7831e-05   Range:   [0,inf),
#[Out]#  Param #3   AMPLITUDE1 =    0.0273161 +/-     0.000413443 ,
#[Out]#  Param #4       SHIFT1 =      127.363 +/-     4.39559e-05 ,
#[Out]#  Param #5       WIDTH1 =   0.00251509 +/-     4.39572e-05   Range:   [0,inf)]
v1 = (sp.specfit.parinfo['SHIFT0'] - 127.36766)/127.36766
v2 = (sp.specfit.parinfo['SHIFT1'] - 127.36766)/127.36766
v1,v2
#[Out]# (-0.00013016989516684947, -3.9638116547612844e-05)
v1 = (sp.specfit.parinfo['SHIFT0'] - 127.36766)/127.36766*constants.c
v2 = (sp.specfit.parinfo['SHIFT1'] - 127.36766)/127.36766*constants.c
v1,v2
from astropy import constants
v1 = (sp.specfit.parinfo['SHIFT0'] - 127.36766)/127.36766*constants.c
v2 = (sp.specfit.parinfo['SHIFT1'] - 127.36766)/127.36766*constants.c
v1,v2
#[Out]# (<Quantity -39023.95282967 m / s>, <Quantity -11883.2083903 m / s>)
v1 = (sp.specfit.parinfo['SHIFT0'] - 127.36766)/127.36766*constants.c.to(u.km/u.s)
v2 = (sp.specfit.parinfo['SHIFT1'] - 127.36766)/127.36766*constants.c.to(u.km/u.s)
v1,v2
#[Out]# (<Quantity -39.02395283 km / s>, <Quantity -11.88320839 km / s>)
v1 = -(sp.specfit.parinfo['SHIFT0'] - 127.36766)/127.36766*constants.c.to(u.km/u.s)
v2 = -(sp.specfit.parinfo['SHIFT1'] - 127.36766)/127.36766*constants.c.to(u.km/u.s)
v1,v2
#[Out]# (<Quantity 39.02395283 km / s>, <Quantity 11.88320839 km / s>)
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.7*u.GHz, xmax=127.8*u.GHz)
#sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
#sp.specfit.parinfo
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.6*u.GHz, xmax=127.8*u.GHz)
#sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
#sp.specfit.parinfo
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.8*u.GHz, xmax=127.9*u.GHz)
#sp.specfit(guesses=[0.03, 127.35, 0.001, 0.003, 127.37, 0.001])
#sp.specfit.parinfo
sp.xarr.convert_to_unit(u.GHz)
sp.plotter(xmin=127.8*u.GHz, xmax=127.9*u.GHz)
sp.specfit(guesses=[0.03, 127.84, 0.001, 0.003, 127.85, 0.001])
#sp.specfit.parinfo
f1 = (sp.specfit.parinfo['SHIFT0'] * (1+ v1/constants.c.to(u.km/u.s)))
f2 = (sp.specfit.parinfo['SHIFT1'] * (1+ v1/constants.c.to(u.km/u.s)))
f1,f2
#[Out]# (<Quantity 127.8568992>, <Quantity 127.86847345>)
f1 = (sp.specfit.parinfo['SHIFT0'] * (1+ v1/constants.c.to(u.km/u.s)))
f2 = (sp.specfit.parinfo['SHIFT1'] * (1+ v2/constants.c.to(u.km/u.s)))
f1,f2
#[Out]# (<Quantity 127.8568992>, <Quantity 127.85689879>)
from astroquery.splatalogue import Splatalogue, utils
utils.minimize_table(Splatalogue.query_lines(f1*(1-5/3e5), f1+(1/3e5)))
from astroquery.splatalogue import Splatalogue, utils
utils.minimize_table(Splatalogue.query_lines(f1*u.GHz*(1-5/3e5), f1+(1/3e5)*u.GHz))
from astroquery.splatalogue import Splatalogue, utils
utils.minimize_table(Splatalogue.query_lines(f1*u.GHz*(1-5/3e5), f1*(1+1/3e5)*u.GHz))
#[Out]# <Table length=41>
#[Out]#      Species          ChemicalName                 QNs                   Freq    log10_Aij    EU_K   
#[Out]#       str18              str17                    str31                float64    float64   float64  
#[Out]# ------------------ ----------------- ------------------------------- ----------- --------- ----------
#[Out]#        AA-n-C4H9CN   n-Butyl cyanide             49(29,20)-48(29,19) 127.8548165  -4.04367  700.45171
#[Out]#        AA-n-C4H9CN   n-Butyl cyanide             49(29,21)-48(29,20) 127.8548165  -4.04367  700.45171
#[Out]#             ClClOO  Chloryl chloride             27(11,16)-27(10,18) 127.8554495  -5.00593   152.8732
#[Out]#           HCCCH2OH Propargyl Alcohol      42(16,26)-41(17,24),vt=1-0 127.8556644  -6.19927  761.77233
#[Out]#           HCCCH2OH Propargyl Alcohol      42(16,26)-41(17,24),vt=1-0 127.8556644  -6.19927  761.77233
#[Out]#             C3H6O2    Hydroxyacetone              17(3,14)-17(2,16)E 127.8558583  -7.57422   57.03782
#[Out]#          CH2CH13CN     Vinyl Cyanide              75(10,65)-76(9,68) 127.8559494  -4.82092 1507.83352
#[Out]#       CHOCHOHCH2OH    Glyceraldehyde             97(26,71)-97(25,73) 127.8560121  -4.78877 1282.25587
#[Out]# NH2CH2CH2OHv25+v27      Aminoethanol               13(1,12)-12(1,11) 127.8561271  -4.09678  633.14529
#[Out]#               C3H8           Propane                             N/A 127.8562559  -8.36421  643.49721
#[Out]#           37ClOv=0 Chlorine monoxide J=7/2-5/2,&Omega;=3/2,F=2-2,l=f  127.856285  -5.55501    10.5242
#[Out]#  CH2OHCOCH2OHv29=1  Dihydroxyacetone             29(16,13)-30(15,16) 127.8563321  -5.64032  386.24425
#[Out]#  CH2OHCOCH2OHv29=1  Dihydroxyacetone             29(16,14)-30(15,15) 127.8563321  -5.64032  386.24425
#[Out]#           37ClOv=0 Chlorine monoxide J=7/2-5/2,&Omega;=3/2,F=2-2,l=e 127.8563487  -5.55501    10.5242
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-2  127.856477       0.0    9.20788
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-2 127.8564858  -5.35999    9.20742
#[Out]#               33SO   Sulfur Monoxide           3(3)-2(2),  F=7/2-7/2  127.856548  -5.50593   25.40698
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-2 127.8565535  -5.36952    9.20731
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-0  127.856571       0.0    9.20788
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-0 127.8565842  -5.01324    9.20742
#[Out]#   CH3CH2CN,v12=1-A     Ethyl Cyanide               21(2,20)-21(1,21) 127.8566351  -5.24138  871.97721
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-0 127.8566481  -5.02277    9.20731
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-1  127.856744       0.0    9.20789
#[Out]#               33SO   Sulfur Monoxide             3(3)-2(2),F=7/2-7/2  127.856767       0.0   25.40843
#[Out]#              CH2NH       Methanimine                   2(0,2)-1(0,1)  127.856769  -4.76741    9.20789
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-1 127.8567698  -4.88289    9.20743
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=3-2  127.856795       0.0    9.20789
#[Out]#              CH2NH       Methanimine                   2(0,2)-1(0,1) 127.8567952  -4.75801    9.20743
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=2-1 127.8568246  -4.89242    9.20732
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=3-2 127.8568268  -4.75801    9.20743
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=3-2 127.8568758  -4.76744    9.20732
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-2  127.856972       0.0     9.2079
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-2 127.8569824  -6.31424    9.20744
#[Out]#        Ga-n-C3H7OH <i>n</i>-Propanol             63(24,40)-62(25,38) 127.8571382  -6.03986 1179.74039
#[Out]#        Ga-n-C3H7OH <i>n</i>-Propanol             63(24,40)-62(25,37) 127.8571382  -5.96756 1179.74039
#[Out]#        Ga-n-C3H7OH <i>n</i>-Propanol             63(24,39)-62(25,38) 127.8571382  -5.96756 1179.74039
#[Out]#        Ga-n-C3H7OH <i>n</i>-Propanol             63(24,39)-62(25,37) 127.8571382  -6.03986 1179.74039
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-1  127.857239       0.0    9.20792
#[Out]#              CH2NH       Methanimine             2(0,2)-1(0,1),F=1-1 127.8572669  -5.13814    9.20745
#[Out]#             CH3C3N   2-Butynenitrile                   31(13)-30(13)   127.85728   -3.6525 1359.83833
#[Out]#             CH3C3N   2-Butynenitrile                   31(13)-30(13)  127.857287  -1.82435 1357.23632
from astroquery.splatalogue import Splatalogue, utils
utils.minimize_table(Splatalogue.query_lines(f1*u.GHz*(1-5/3e5), f1*(1+5/3e5)*u.GHz))
#[Out]# <Table length=56>
#[Out]#      Species             ChemicalName                        QNs                      Freq    log10_Aij    EU_K   
#[Out]#       str18                 str24                           str37                   float64    float64   float64  
#[Out]# ------------------ ------------------------ ------------------------------------- ----------- --------- ----------
#[Out]#        AA-n-C4H9CN          n-Butyl cyanide                   49(29,20)-48(29,19) 127.8548165  -4.04367  700.45171
#[Out]#        AA-n-C4H9CN          n-Butyl cyanide                   49(29,21)-48(29,20) 127.8548165  -4.04367  700.45171
#[Out]#             ClClOO         Chloryl chloride                   27(11,16)-27(10,18) 127.8554495  -5.00593   152.8732
#[Out]#           HCCCH2OH        Propargyl Alcohol            42(16,26)-41(17,24),vt=1-0 127.8556644  -6.19927  761.77233
#[Out]#           HCCCH2OH        Propargyl Alcohol            42(16,26)-41(17,24),vt=1-0 127.8556644  -6.19927  761.77233
#[Out]#             C3H6O2           Hydroxyacetone                    17(3,14)-17(2,16)E 127.8558583  -7.57422   57.03782
#[Out]#          CH2CH13CN            Vinyl Cyanide                    75(10,65)-76(9,68) 127.8559494  -4.82092 1507.83352
#[Out]#       CHOCHOHCH2OH           Glyceraldehyde                   97(26,71)-97(25,73) 127.8560121  -4.78877 1282.25587
#[Out]# NH2CH2CH2OHv25+v27             Aminoethanol                     13(1,12)-12(1,11) 127.8561271  -4.09678  633.14529
#[Out]#               C3H8                  Propane                                   N/A 127.8562559  -8.36421  643.49721
#[Out]#           37ClOv=0        Chlorine monoxide       J=7/2-5/2,&Omega;=3/2,F=2-2,l=f  127.856285  -5.55501    10.5242
#[Out]#  CH2OHCOCH2OHv29=1         Dihydroxyacetone                   29(16,13)-30(15,16) 127.8563321  -5.64032  386.24425
#[Out]#  CH2OHCOCH2OHv29=1         Dihydroxyacetone                   29(16,14)-30(15,15) 127.8563321  -5.64032  386.24425
#[Out]#           37ClOv=0        Chlorine monoxide       J=7/2-5/2,&Omega;=3/2,F=2-2,l=e 127.8563487  -5.55501    10.5242
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-2  127.856477       0.0    9.20788
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-2 127.8564858  -5.35999    9.20742
#[Out]#               33SO          Sulfur Monoxide                 3(3)-2(2),  F=7/2-7/2  127.856548  -5.50593   25.40698
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-2 127.8565535  -5.36952    9.20731
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-0  127.856571       0.0    9.20788
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-0 127.8565842  -5.01324    9.20742
#[Out]#   CH3CH2CN,v12=1-A            Ethyl Cyanide                     21(2,20)-21(1,21) 127.8566351  -5.24138  871.97721
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-0 127.8566481  -5.02277    9.20731
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-1  127.856744       0.0    9.20789
#[Out]#               33SO          Sulfur Monoxide                   3(3)-2(2),F=7/2-7/2  127.856767       0.0   25.40843
#[Out]#              CH2NH              Methanimine                         2(0,2)-1(0,1)  127.856769  -4.76741    9.20789
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-1 127.8567698  -4.88289    9.20743
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=3-2  127.856795       0.0    9.20789
#[Out]#              CH2NH              Methanimine                         2(0,2)-1(0,1) 127.8567952  -4.75801    9.20743
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=2-1 127.8568246  -4.89242    9.20732
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=3-2 127.8568268  -4.75801    9.20743
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=3-2 127.8568758  -4.76744    9.20732
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-2  127.856972       0.0     9.2079
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-2 127.8569824  -6.31424    9.20744
#[Out]#        Ga-n-C3H7OH        <i>n</i>-Propanol                   63(24,40)-62(25,38) 127.8571382  -6.03986 1179.74039
#[Out]#        Ga-n-C3H7OH        <i>n</i>-Propanol                   63(24,40)-62(25,37) 127.8571382  -5.96756 1179.74039
#[Out]#        Ga-n-C3H7OH        <i>n</i>-Propanol                   63(24,39)-62(25,38) 127.8571382  -5.96756 1179.74039
#[Out]#        Ga-n-C3H7OH        <i>n</i>-Propanol                   63(24,39)-62(25,37) 127.8571382  -6.03986 1179.74039
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-1  127.857239       0.0    9.20792
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-1 127.8572669  -5.13814    9.20745
#[Out]#             CH3C3N          2-Butynenitrile                         31(13)-30(13)   127.85728   -3.6525 1359.83833
#[Out]#             CH3C3N          2-Butynenitrile                         31(13)-30(13)  127.857287  -1.82435 1357.23632
#[Out]#              CH2NH              Methanimine                   2(0,2)-1(0,1),F=1-1  127.857326  -5.14767    9.20734
#[Out]#   NH2CH2CH2OHv25=1             Aminoethanol                      16(9,8)-17(7,11) 127.8575339  -7.84786  541.73578
#[Out]#             l-HC6N Cyanobutadiynylmethylene                       N=76-75,J=75-74 127.8576142  -3.49423  236.72129
#[Out]#             O79BrO          Bromine Dioxide 23(5,18)-23(4,19),J=47/2-47/2,F=23-23 127.8576844  -4.13085  221.76184
#[Out]#             O79BrO          Bromine Dioxide 23(5,18)-23(4,19),J=47/2-47/2,F=23-24  127.857809  -6.58245  221.76184
#[Out]#             O79BrO          Bromine Dioxide 23(5,18)-23(4,19),J=47/2-47/2,F=22-22 127.8580314  -4.12927  221.76013
#[Out]#               33SO          Sulfur Monoxide                   3(3)-2(2),F=3/2-5/2  127.858309       0.0    25.4085
#[Out]#        AG-n-C4H9CN          n-Butyl cyanide                   44(32,12)-43(32,11) 127.8583702  -4.11031 1088.19472
#[Out]#        AG-n-C4H9CN          n-Butyl cyanide                   44(32,13)-43(32,12) 127.8583702  -4.11031 1088.19472
#[Out]#  aG'g-CH3CHOHCH2OH    1,2-propanediol, aG'g                     61(4,57)-60(5,55) 127.8585279  -7.16481    561.678
#[Out]#  aG'g-CH3CHOHCH2OH    1,2-propanediol, aG'g                     61(5,57)-60(6,55) 127.8585279  -7.16481    561.678
#[Out]#             O79BrO          Bromine Dioxide 23(5,18)-23(4,19),J=47/2-47/2,F=24-23 127.8588905  -6.58935  221.76189
#[Out]#               33SO          Sulfur Monoxide                 3(3)-2(2),  F=3/2-5/2 127.8590128  -6.35969   25.40825
#[Out]#             O79BrO          Bromine Dioxide 23(5,18)-23(4,19),J=47/2-47/2,F=24-24 127.8590152  -4.13075   221.7619
#[Out]#           H2C(CN)2            Malononitrile                     31(9,23)-32(8,24)  127.859029  -4.63728  202.72904
import sqlite3
conn = sqlite3.connect('/orange/adamginsburg/software/XCLASS-Interface/Database/cdms_sqlite.db')
conn.cursor
#[Out]# <function Connection.cursor>
conn.cursor()
#[Out]# <sqlite3.Cursor at 0x2b21f12e1340>
cursor = conn.cursor()
cursor.execute("SELECT PF_ID, PF_Name, PF_SpeciesID, "
                           "PF_VamdcSpeciesID, PF_Recommendation, "
                           "PF_Status, PF_Createdate, "
                           "PF_Checkdate FROM Partitionfunctions "
                           "WHERE PF_Status = ?", (status,))
cursor = conn.cursor()
            cursor.execute("SELECT PF_ID, PF_Name, PF_SpeciesID, PF_VamdcSpeciesID, \
                            PF_Recommendation, PF_Status, PF_Createdate, \
                            PF_Checkdate FROM Partitionfunctions")
cursor = conn.cursor()
cursor.execute("SELECT PF_ID, PF_Name, PF_SpeciesID, PF_VamdcSpeciesID, \
                PF_Recommendation, PF_Status, PF_Createdate, \
                PF_Checkdate FROM Partitionfunctions")
cursor = conn.cursor()
cursor.execute("SELECT *")
cursor = conn.cursor()
cursor.execute("SELECT * FROM *")
cursor = conn.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(cursor.fetchall())
cursor = conn.cursor()
cursor.execute("SELECT * FROM Transitions")
print(cursor.fetchall())
