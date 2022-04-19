########################################################
# Started Logging At: 2022-03-21 17:12:54
########################################################
########################################################
# # Started Logging At: 2022-03-21 17:12:55
########################################################
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'

# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
from astropy import coordinates, units as u
from astropy import wcs
from astropy.io import fits
crd = coordinates.SkyCoord("17:46:10.80 -28:42:13.3", frame='fk5', unit=(u.h, u.deg))
import warnings
warnings.filterwarnings('ignore')
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
tempmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/temperature_map_februrary.fits')
colmap = fits.open('/blue/adamginsburg/abulatek/brick/first_results/temperature_map/col_density_map_februrary.fits')
ww = wcs.WCS(temp[0].header)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        restfrq = ch3cntbl[np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])]['Freq-GHz(rest frame,redshifted)']
        print(fn, restfrq)
        #ccube = ch3cnube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, 
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = ch3cntbl[qns][0]
        print(fn, qns, restfrq)
        ccube = ch3cnube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 250*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = ch3cntbl[msk][0]
        print(fn, qns, restfrq)
        ccube = ch3cnube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 250*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = ch3cntbl[msk][0]
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 250*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = ch3cntbl[msk]['Resolved QNs'][0]
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 250*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 250*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
        ccube.write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
        ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.fits', overwrite=True)
get_ipython().run_line_magic('ls', '*pbcor* -d')
get_ipython().run_line_magic('ls', '*pbcor* $results/-d')
get_ipython().run_line_magic('ls', '$results/*pbcor* -d')
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
cubefns = glob.glob(f"{results}/*.image.pbcor")
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
    if any('(0)' in x for x in ch3cntbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cntbl['Resolved QNs']])
        restfrq = ch3cntbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cntbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        print(fn, qns, restfrq)
        ccube = ch3cncube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
        ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cn_cutouts/CH3CN_{qns}.pbcor.fits', overwrite=True)
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
