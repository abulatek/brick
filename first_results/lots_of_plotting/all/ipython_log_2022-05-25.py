########################################################
# Started Logging At: 2022-05-25 15:09:06
########################################################
########################################################
# # Started Logging At: 2022-05-25 15:09:07
########################################################
########################################################
# Started Logging At: 2022-05-25 15:33:40
########################################################
########################################################
# # Started Logging At: 2022-05-25 15:33:40
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
ind = 25
from spectral_cube import SpectralCube

cube = SpectralCube.read(results+'source_ab_'+freq_spws[ind]+'_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image', use_dask = True).to(u.K)
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
# Extract spectrum for brightest pixel in central core region
import regions
reg = regions.Regions.read('/blue/adamginsburg/abulatek/brick/first_results/region_temperatures/centralcoreellipse.reg')

smoothed_cube_reg = smoothed_cube.subcube_from_regions(reg)
m0 = smoothed_cube_reg.moment0()

import numpy as np

brightest_pixel = np.unravel_index(np.nanargmax(m0), m0.shape)
spectrum = smoothed_cube_reg[:, brightest_pixel[0], brightest_pixel[1]]
import matplotlib.pyplot as plt

SMALL_SIZE = 10
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
fig = plt.figure(figsize = (30, 8))

plt.plot(smoothed_cube_reg.spectral_axis.to(u.GHz), np.array(spectrum), linestyle = '-', color = 'k', linewidth = 1, drawstyle = 'steps-mid')
plt.xlabel(f"Frequency [{smoothed_cube_reg.spectral_axis.to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
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
source_size =  0.33 # total guess, but based on Dan's data? also, not used
tkin = 153.1 # 7.9842254347e+01
Ntot = 1.5e+14 # 10**N_tot # 1.5061155735e+14 # 1.0061155735e+16
vwidth_fwhm = 4.2101396644e+00
vwidth = vwidth_fwhm / np.sqrt(8*np.log(2))
vcen = 38.0
tbg = 0
import sys
import os
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

# depending on parameter "Inter_Flag" define beam size (in arcsec)
# (Inter_Flag = True) or size of telescope (in m) (Inter_Flag = False)
TelescopeSize = 1000.0 # meters

# define beam minor axis length (in arsec)
BMIN = smoothed_cube_reg.beam.minor.to(u.arcsec).value

# define beam major axis length (in arsec)
BMAJ = smoothed_cube_reg.beam.major.to(u.arcsec).value

# define beam position angle (in degree)
BPA = smoothed_cube_reg.beam.pa.to(u.deg).value

# interferrometric data?
Inter_Flag = True # SCREAM (i.e., maybe check this?)

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

# set_xclass_parameters() # The molfit file is already created, I don't want to mess with it; to run this, need to uncomment function definition a few cells ago
modeldata, log, TransEnergies, IntOpt, JobDir = myxclass_call(FreqMin=FreqMin, FreqMax=FreqMax)
fig = plt.figure(figsize = (30, 8))

plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], linestyle = '-', color = 'r', linewidth = 1, drawstyle = 'steps-mid')
plt.xlabel(f"Frequency [{(modeldata[:,0]*u.MHz).to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
plt.show()
# Preliminary continuum subtraction
spectrum_contsub = np.array(spectrum) - np.median(np.array(spectrum))
fig = plt.figure(figsize = (30, 8))

plt.plot(cube.spectral_axis.to(u.GHz), spectrum_contsub, linestyle = '-', color = 'k', linewidth = 0.75, drawstyle = 'steps-mid', 
         label = "Data")
plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], linestyle = '--', color = 'r', linewidth = 1, drawstyle = 'steps-mid', 
         label = "Model")
plt.legend(loc='best')
plt.xlabel(f"Frequency [{smoothed_cube_reg.spectral_axis.to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
plt.title(f"{freq_spws[ind]} data and model comparison, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92
plt.show()
import glob
MolfitsFileNames = glob.glob(LocalPath + 'individuals/*.molfit')
modeldata_all = []
transenergies_all = []
for MolfitsFileName in MolfitsFileNames:
    modeldata, log, TransEnergies, IntOpt, JobDir = myxclass_call(FreqMin=FreqMin, FreqMax=FreqMax)
    modeldata_all.append(modeldata)
    transenergies_all.append(TransEnergies)
fig = plt.figure(figsize = (30, 8))

plt.plot(cube.spectral_axis.to(u.GHz), spectrum_contsub, linestyle = '-', color = 'k', linewidth = 1, drawstyle = 'steps-mid', 
        label = "Data")
for modeldata, transenergies in zip(modeldata_all, transenergies_all):
    if len(modeldata) != 0:
        plt.plot((modeldata[:,0]*u.MHz).to(u.GHz), modeldata[:,1], linestyle = '--', linewidth = 1, drawstyle = 'steps-mid', 
                 label = transenergies[1][-1])
    else:
        print("***** Warning! One of the models is a blank list! *****")
plt.legend(loc='best')
plt.xlabel(f"Frequency [{smoothed_cube_reg.spectral_axis.to(u.GHz).unit}]")
plt.ylabel(f"Brightness temperature [{spectrum.unit}]")
plt.title(f"{freq_spws[ind]} data and model comparison, $T$ = {tkin:.1f} K, $\log_{{10}}(N_{{tot}})$ = {np.log10(Ntot):.2f}") # y = 0.92
plt.show()
smoothed_cube_reg.beam
#[Out]# Beam: BMAJ=1.5765712380621564 arcsec BMIN=1.2710207774915303 arcsec BPA=94.07396509688701 deg
########################################################
# Started Logging At: 2022-05-25 18:16:06
########################################################
########################################################
# # Started Logging At: 2022-05-25 18:16:07
########################################################
