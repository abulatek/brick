########################################################
# Started Logging At: 2022-09-14 16:13:14
########################################################
########################################################
# # Started Logging At: 2022-09-14 16:13:15
########################################################
########################################################
# Started Logging At: 2022-09-14 16:14:23
########################################################
########################################################
# # Started Logging At: 2022-09-14 16:14:24
########################################################
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
import numpy as np
# %matplotlib inline
from spectral_cube import SpectralCube
from astropy import units as u
results = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
cubefns = glob.glob(f"{results}/*.image.pbcor")
from astropy import coordinates, units as u
from astropy import wcs
from astropy.io import fits
crd = coordinates.SkyCoord("17:46:10.80 -28:42:13.3", frame='fk5', unit=(u.h, u.deg))
cubes = []

for fn in cubefns:

    ch3cchcube = SpectralCube.read(fn, format='casa_image')
    
    ch3cchtbl = ch3cchcube.find_lines(chemical_name='CH3CCH', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    #ch3cchtbl = ch3cchtbl[ch3cchtbl['Quantum Number Code'] == 202]
    if len(ch3cchtbl) > 0:
        cubes.append(ch3cchcube)
        print(fn)
        print(ch3cchtbl)
    if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
        restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        for rf, qn in zip(restfrq, qns):
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            #ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cch_cutouts/CH3CCH_{qn}.pbcor.fits', overwrite=True)
import glob
cubefns = glob.glob(f"{results}/source_ab_*_clean_2sigma_n50000_masked_3sigma_pbmask0p18.image")
#cubefns = glob.glob(f"{results}/*.image.pbcor")
from astropy import coordinates, units as u
from astropy import wcs
from astropy.io import fits
crd = coordinates.SkyCoord("17:46:10.80 -28:42:13.3", frame='fk5', unit=(u.h, u.deg))
import warnings
warnings.filterwarnings('ignore')
cubes = []

for fn in cubefns:

    ch3cchcube = SpectralCube.read(fn, format='casa_image')
    
    ch3cchtbl = ch3cchcube.find_lines(chemical_name='CH3CCH', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    #ch3cchtbl = ch3cchtbl[ch3cchtbl['Quantum Number Code'] == 202]
    if len(ch3cchtbl) > 0:
        cubes.append(ch3cchcube)
        print(fn)
        print(ch3cchtbl)
    if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
        restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        for rf, qn in zip(restfrq, qns):
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            #ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cch_cutouts/CH3CCH_{qn}.pbcor.fits', overwrite=True)
cubes = []

for fn in cubefns:

    ch3cchcube = SpectralCube.read(fn, format='casa_image')
    print(fn)
    
    ch3cchtbl = ch3cchcube.find_lines(chemical_name='^CH3CCH', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    #ch3cchtbl = ch3cchtbl[ch3cchtbl['Quantum Number Code'] == 202]
    if len(ch3cchtbl) > 0:
        cubes.append(ch3cchcube)
        print(fn)
        print(ch3cchtbl)
    if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
        restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        for rf, qn in zip(restfrq, qns):
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            #ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cch_cutouts/CH3CCH_{qn}.pbcor.fits', overwrite=True)
cubes = []

for fn in cubefns:

    ch3cchcube = SpectralCube.read(fn, format='casa_image')
    print(fn)
    
    ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)
    #ch3cchtbl = ch3cchtbl[ch3cchtbl['Quantum Number Code'] == 202]
    if len(ch3cchtbl) > 0:
        cubes.append(ch3cchcube)
        print(fn)
        print(ch3cchtbl)
    if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
        msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
        restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
        qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
        for rf, qn in zip(restfrq, qns):
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            #ccube.to(u.K).write(f'/orange/adamginsburg/brick_alma_linesurvey/ch3cch_cutouts/CH3CCH_{qn}.pbcor.fits', overwrite=True)
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)

        x,y = map(int, ch3cchcube.wcs.celestial.world_to_pixel(crd))

        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        

        data_sp = ch3cchcube[:, y, x]
        data_sp_K = data_sp.value * ch3cchcube.jtok_factors()
        ax.plot(ch3cchcube.spectral_axis.to(u.GHz), data_sp_K)
        ax.set_title(ch3cchtbl['Resolved QNs'][0])
print("We did it!")
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)

        x,y = map(int, ch3cchcube.wcs.celestial.world_to_pixel(crd))

        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis.to(u.GHz), data_sp_K)
                ax.set_title(ccube['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis.to(u.GHz), data_sp_K)
                ax.set_title(ccube['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
                print(ccube)

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis.to(u.GHz), data_sp_K)
                ax.set_title(ccube['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
                if ccube.shape[0] == 1:
                    continue # ?!?!?!? why

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis.to(u.GHz), data_sp_K)
                ax.set_title(ccube['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
                if ccube.shape[0] == 1:
                    continue # ?!?!?!? why

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis, data_sp_K)
                ax.set_title(ccube['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
                if ccube.shape[0] == 1:
                    continue # ?!?!?!? why

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis, data_sp_K)
                ax.set_title(ch3cchtbl['Resolved QNs'][0])
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:
        ax = pl.subplot(3,3,ind+1)


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            for rf, qn in zip(restfrq, qns):
                print(fn, qn, rf)
                ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
                x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
                if ccube.shape[0] == 1:
                    continue # ?!?!?!? why

                data_sp = ccube[:, y, x]
                data_sp_K = data_sp.value * ccube.jtok_factors()
                ax.plot(ccube.spectral_axis, data_sp_K)
                ax.set_title(qn)
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            print(fn, qn, rf)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
cubes
#[Out]# [DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 271307763969.477 Hz:273182245597.862 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 135285431368.900 Hz:137159997643.322 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 273010861630.481 Hz:274885343259.334 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 255308020805.109 Hz:257182502433.963 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 152395310730.523 Hz:154269877084.026 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(1920, 512, 512) and unit=Jy / beam and chunk size (80, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   1920  type_s: FREQ      unit_s: Hz     range: 101549443159.977 Hz:103423216948.801 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 255308020805.109 Hz:257182502433.963 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(3840, 512, 512) and unit=Jy / beam and chunk size (128, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   3840  type_s: FREQ      unit_s: Hz     range: 255308020805.109 Hz:257182502433.963 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(1920, 512, 512) and unit=Jy / beam and chunk size (80, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   1920  type_s: FREQ      unit_s: Hz     range: 101915892799.071 Hz:103789669232.376 Hz,
#[Out]#  DaskVaryingResolutionSpectralCube with shape=(1920, 512, 512) and unit=Jy / beam and chunk size (80, 128, 512):
#[Out]#   n_x:    512  type_x: RA---SIN  unit_x: deg    range:   266.528130 deg:  266.560501 deg
#[Out]#   n_y:    512  type_y: DEC--SIN  unit_y: deg    range:   -28.719152 deg:  -28.690763 deg
#[Out]#   n_s:   1920  type_s: FREQ      unit_s: Hz     range: 101549438071.648 Hz:103423209520.988 Hz]
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0]).replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.figure()
#[Out]# <Figure size 432x288 with 0 Axes>
pl.show()
get_ipython().run_line_magic('matplotlib', 'inline')
pl.show()
pl.plot([0,1])
#[Out]# [<matplotlib.lines.Line2D at 0x2b35037bfe20>]
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.show())
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.show()
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=rf*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.show()
pl.figure(figsize=(12,12))
for ind, ch3cchcube in enumerate(cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 350*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
pl.show()
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            print("theory",restfrq)
            restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
            print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 200*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
            ind += 1
pl.show()
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            print("theory",restfrq)
            if restfrq.mask:
                restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
                print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 200*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
            ind += 1
pl.show()
restfrq
#[Out]# array([273.4199209])
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 9:


        # TODO: determine this from data
        v_cen = 10*u.km/u.s
        v_disp = 1.5*u.km/u.s
        
        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            print("theory",restfrq)
            if not np.all(np.isfinite(restfrq)):
                restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
                print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-50*u.km/u.s, 200*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
            ind += 1
pl.show()
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 16:

        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            print("theory",restfrq)
            if not np.all(np.isfinite(restfrq)):
                restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
                print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-30*u.km/u.s, 150*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(4,4,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K)
            ax.set_title(str(qns))
            ind += 1
pl.show()
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 16:

        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            print("theory",restfrq)
            if not np.all(np.isfinite(restfrq)):
                restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
                print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-30*u.km/u.s, 150*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(4,4,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K, drawstyle='steps-mid', color='k')
            ax.set_title(str(qns))
            ind += 1
pl.show()
pl.tight_layout()
pl.figure(figsize=(12,12))
ind = 0
for ch3cchcube in (cubes):
    if ind < 9:

        ch3cchtbl = ch3cchcube.find_lines(chemical_name=' CH3CCH v = 0', line_lists=['JPL'], 
                                    show_upper_degeneracy=True, show_qn_code=True)        
        
        if any('(0)' in x for x in ch3cchtbl['Resolved QNs']):
            msk = np.array(['(0)' in x for x in ch3cchtbl['Resolved QNs']])
            restfrq = ch3cchtbl[msk]['Freq-GHz(rest frame,redshifted)'].value
            #print("theory",restfrq)
            if not np.all(np.isfinite(restfrq)):
                restfrq = ch3cchtbl[msk]['Meas Freq-GHz(rest frame,redshifted)'].value
                #print("meas", restfrq)
            qns = str(ch3cchtbl[msk]['Resolved QNs'][0])#.replace("(","").replace(")","")
            #print(restfrq, qns)
            ccube = ch3cchcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.GHz, velocity_convention='radio').spectral_slab(-30*u.km/u.s, 150*u.km/u.s)
            x,y = map(int, ccube.wcs.celestial.world_to_pixel(crd))
            if ccube.shape[0] == 1:
                print("NOPE")
                continue # ?!?!?!? why

            ax = pl.subplot(3,3,ind+1)
            data_sp = ccube[:, y, x]
            data_sp_K = data_sp.value * ccube.jtok_factors()
            ax.plot(ccube.spectral_axis, data_sp_K, drawstyle='steps-mid', color='k')
            ax.set_title(str(qns))
            ind += 1
pl.show()
pl.tight_layout()
