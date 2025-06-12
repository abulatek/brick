from spectral_cube import SpectralCube
from astropy.io import fits
from reproject import reproject_interp
from astroquery.splatalogue.utils import minimize_table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import pyspeckit
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.lte_molecule import nupper_of_kkms
from astropy import coordinates, constants, units as u
from astropy.stats import mad_std

def fetch_cubes(cubefns, catalog, mol_name_pyspeckit=None, mol_tag_pyspeckit=None, parse_loc=False, 
                ret_tbl=False):
    # Get only the cubes that have our line of interest in them
    cubes = []
    tbls = []
    for fn in cubefns:
        molcube = SpectralCube.read(fn, format='fits', use_dask=False) # My version of dask breaks this so bad
        # print(molcube)
        try:
            if ret_tbl:
                freqs, aij, deg, EU, partfunc, tbl = get_molecular_parameters(molecule_name=mol_name_pyspeckit,
                                                                              molecule_tag=mol_tag_pyspeckit,
                                                                              catalog=catalog, 
                                                                              parse_name_locally=parse_loc, 
                                                                              return_table=ret_tbl,
                                                                              fmin=molcube.spectral_axis.min(), 
                                                                              fmax=molcube.spectral_axis.max())
            else: 
                freqs, aij, deg, EU, partfunc = get_molecular_parameters(molecule_name=mol_name_pyspeckit,
                                                                         molecule_tag=mol_tag_pyspeckit, 
                                                                         catalog=catalog, 
                                                                         parse_name_locally=parse_loc, 
                                                                         return_table=ret_tbl,
                                                                         fmin=molcube.spectral_axis.min(), 
                                                                         fmax=molcube.spectral_axis.max())
        except TypeError as te:
            # print(f"ERROR: {te}")
            continue  
        except:
            continue
        if len(freqs) > 0:
            cubes.append(molcube)
            # print(molcube.beam)
            if ret_tbl:
                tbls.append(tbl)
    # Reorder selected cubes from lowest to highest frequency
    lowestfreqs = [np.min(molcube.spectral_axis).value for molcube in cubes]
    correct_order = np.argsort(lowestfreqs)
    cubes_inorder = [cubes[i] for i in correct_order]
    # Return cubes (and tables if requested)
    if ret_tbl: 
        # Reorder tables from lowest to highest frequency
        lowestfreqs = [np.min(tbl['FREQ']) for tbl in tbls]
        correct_order = np.argsort(lowestfreqs)
        tbls_inorder = [tbls[i] for i in correct_order]
        print([tbl for tbl in tbls_inorder])
        print("\n*** NOTE: Make sure you are capturing the cubes and tables into different variables!")
        return cubes_inorder, tbls_inorder
    else:
        return cubes_inorder

def get_cubes_from_mask(region_fn, region_ind, cubes, plot_region=False, return_mask=False):
    # Set region with diffuse emission
    diffuse_regions = fits.open(region_fn)
    diffuse_region_mask = diffuse_regions[0].data == region_ind
    # Plot which region you are using (optional)
    if plot_region == True:
        plt.imshow(diffuse_region_mask, origin="lower")
    cubes_masked = []
    for molcube in cubes:
        array, footprint = reproject_interp((diffuse_region_mask, diffuse_regions[0].header), 
                                             molcube.wcs.celestial, shape_out=(molcube.shape[1],
                                                                               molcube.shape[2]))
        mask = array == 1
        # Get subcube from the mask
        molcube_masked = molcube.with_mask(mask)
        cubes_masked.append(molcube_masked)
    if not return_mask:
        return cubes_masked
    if return_mask:
        return cubes_masked, mask

def model_and_plot(cubes, temp, N_tot, v_cen, v_disp, catalog, fig_width, fig_height, nrows, ncols, 
                   mol_name_pyspeckit=None, mol_tag_pyspeckit=None, parse_loc=True, ret_tbl=True, 
                   line_by_line=False, name_for_plot=None, print_diag=False, extr_type=None, crd=None, 
                   just_data=False, EU_cutoff_K=100, aij_cutoff=-5, LGINT_cutoff=-5, calc_N_uppers=False, 
                   show_2_sigma=False, return_freqs=False):
    # Plot all lines of the molecule and overplot model spectrum, given an input T, N_tot, and velocity params
    # Assumes input cubes are continuum-subtracted
    fig = plt.figure(figsize = (fig_width, fig_height))
    # Formatting for plot title
    if extr_type == "coord":
        extr_type_plot = "core"
    elif extr_type == "reg":
        extr_type_plot = "frown"
    plt.suptitle(f"{name_for_plot}, {extr_type_plot}, $T_{{rot}} =$ {temp:0.1f} K, $N_{{tot}} =$ {N_tot:0.5g} cm$^{{-2}}$, $v_{{cen}} =$ {v_cen.value} km s$^{{-1}}$, $v_{{disp}} =$ {v_disp.value} km s$^{{-1}}$", fontsize = 16)
    ind=0 # Need this for subplots
    if calc_N_uppers:
        EUs = []
        log_N_upper_gs = []
    if return_freqs:
        mom0_freqs = []
    # Loop through cubes and extract lines
    for molcube in cubes:
        try:
            if ret_tbl:
                freqs, aij, deg, EU, partfunc, tbl = get_molecular_parameters(molecule_name=mol_name_pyspeckit, 
                                                                              molecule_tag=mol_tag_pyspeckit, 
                                                                              catalog=catalog, 
                                                                              parse_name_locally=parse_loc, 
                                                                              return_table=ret_tbl, 
                                                                              fmin=molcube.spectral_axis.min(), 
                                                                              fmax=molcube.spectral_axis.max())
            else: 
                freqs, aij, deg, EU, partfunc = get_molecular_parameters(molecule_name=mol_name_pyspeckit, 
                                                                         molecule_tag=mol_tag_pyspeckit, 
                                                                         catalog=catalog, 
                                                                         parse_name_locally=parse_loc, 
                                                                         return_table=ret_tbl, 
                                                                         fmin=molcube.spectral_axis.min(), 
                                                                         fmax=molcube.spectral_axis.max())
        except TypeError as te:
            # print(f"ERROR: {te}")
            continue  
        except:
            continue
        if print_diag: print("Got molecular parameters")
        # Do line extraction on a spw-by-spw basis
        if not line_by_line:
            # Do a temporary extraction to get frequency boundaries
            if extr_type == "coord":
                # Extract spectrum from provided coordinate
                x,y = map(int, molcube.wcs.celestial.world_to_pixel(crd))
                molcube_temp = molcube[:, y, x]
            elif extr_type == "reg":
                # Extract spectrum from the subcube made from the mask
                molcube_temp = molcube.mean(axis=(1,2))
            else:
                raise ValueError("Invalid extraction type specification")
            # Make a temporary model spectrum
            mod_temp = lte_molecule.generate_model(molcube_temp.spectral_axis, v_cen, v_disp, temp, N_tot, 
                                                   freqs, aij, deg, EU, partfunc)
            # If there is signal in the model spectrum, get the boundaries around the lines
            if sum(mod_temp > 0.01*mod_temp.max()) > 0:       
                ax = plt.subplot(nrows, ncols, ind+1)
                minfreq = molcube.spectral_axis[mod_temp > 0.01*mod_temp.max()].min().to(u.GHz) - 0.01*u.GHz
                maxfreq = molcube.spectral_axis[mod_temp > 0.01*mod_temp.max()].max().to(u.GHz) + 0.01*u.GHz
                if print_diag: print("Calculated spectral line boundaries")
            else:
                continue
            # Extract the spectrum
            molcube_slice = molcube.spectral_slab(minfreq, maxfreq)
            if extr_type == "coord":
                # Extract spectrum from provided coordinate
                x,y = map(int, molcube.wcs.celestial.world_to_pixel(crd))
                data_sp = molcube_slice[:, y, x]
            elif extr_type == "reg":
                # Extract spectrum from the subcube made from the mask
                data_sp = molcube_slice.mean(axis=(1,2))
            else:
                raise ValueError("Invalid extraction type specification")
            # Continuum-subtract the spectrum
            # data_sp_contsub = data_sp - np.median(data_sp)
            data_sp_contsub = data_sp # Assumes contsub is already done
            if print_diag: print("extracted continuum-subtracted spectrum")
            # Make the model spectrum
            mod = lte_molecule.generate_model(molcube_slice.spectral_axis, v_cen, v_disp, temp, N_tot, 
                                              freqs, aij, deg, EU, partfunc)
            if print_diag: print("made model")
            # Plot the data and model spectra
            ax.plot(molcube_slice.spectral_axis.to(u.GHz), data_sp_contsub, color='k')
            if not just_data:
                ax.plot(molcube_slice.spectral_axis.to(u.GHz), mod)
            ax.set_xlabel(f"Frequency [{(molcube_slice.spectral_axis.to(u.GHz)).unit}]")
            ax.set_ylabel(f"Brightness temperature [{data_sp_contsub.unit}]")
            EU_K = ((EU*u.erg)/constants.k_B).decompose()
            ax.text(0.05, 0.9, f"$E_U = ${EU_K[-1]:0.1f}", transform=ax.transAxes)
            if ret_tbl:
                ax.text(0.05, 0.8, f"{tbl['Ju'][-1]}({tbl['Ku'][-1]}) – {tbl['Jl'][-1]}({tbl['Kl'][-1]})", 
                        transform=ax.transAxes)
            # Show the local 2-sigma mark on each axis if requested (for upper-limit measurements)
            if show_2_sigma:
                # Get noise level for spectrum
                noise_level = mad_std(data_sp_contsub.data)*u.K
                ax.hlines(y=2*noise_level.value, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], linewidth=1, 
                          color='r')
            ind+=1
        # Do line extraction on a line-by-line basis
        elif line_by_line:
            for freq, this_EU, this_aij, this_deg, this_tbl in zip(freqs, EU, aij, deg, tbl):
                this_EU_K = ((this_EU*u.erg)/constants.k_B).decompose()
                # print(freq, this_EU_K, this_aij)
                this_LGINT = this_tbl['LGINT']
                if this_EU_K < EU_cutoff_K*u.K and this_LGINT > LGINT_cutoff: # this_aij > aij_cutoff 
                    ax = plt.subplot(nrows, ncols, ind+1)
                    # Set center at vcen, then do bandwidth of 50 km/s
                    freq_to_vel = u.doppler_radio(freq)
                    minfreq = (v_cen - 50*u.km/u.s).to(u.GHz, equivalencies=freq_to_vel)
                    maxfreq = (v_cen + 50*u.km/u.s).to(u.GHz, equivalencies=freq_to_vel)
                    if print_diag: print("calculated spectral line boundaries")
                    # Extract the spectrum
                    molcube_slice = molcube.spectral_slab(minfreq, maxfreq)
                    if extr_type == "coord":
                        # Extract spectrum from provided coordinate
                        x,y = map(int, molcube.wcs.celestial.world_to_pixel(crd))
                        data_sp = molcube_slice[:, y, x]
                    elif extr_type == "reg":
                        # Extract spectrum from the subcube made from the mask
                        data_sp = molcube_slice.mean(axis=(1,2))
                    else:
                        raise ValueError("Invalid extraction type specification")
                    # Continuum-subtract the spectrum
                    # data_sp_contsub = data_sp - np.median(data_sp)
                    data_sp_contsub = data_sp # Assumes contsub is already done
                    if print_diag: print("extracted continuum-subtracted spectrum")
                    # Make the model spectrum
                    mod = lte_molecule.generate_model(molcube_slice.spectral_axis, v_cen, v_disp, temp, N_tot,
                                                      freqs, aij, deg, EU, partfunc)
                    if print_diag: print("made model")
                    # Plot the data and model spectra
                    ax.plot(molcube_slice.spectral_axis.to(u.GHz), data_sp_contsub, color='k')
                    if not just_data:
                        ax.plot(molcube_slice.spectral_axis.to(u.GHz), mod)
                    ax.ticklabel_format(useOffset=False, style='plain')
                    fig.supxlabel(f"Frequency [{(molcube_slice.spectral_axis.to(u.GHz)).unit}]", y=0.05)
                    fig.supylabel(f"Brightness temperature [{data_sp_contsub.unit}]", x=0.01)
                    if print_diag: print(this_EU_K)
                    # Add labels
                    handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 3
                    labels = [f"$E_U =${this_EU_K:0.2f}", f"log$_{{10}}(A_{{i,j}}) =${str(round(this_aij, 2)).replace('-','−')}", f"LGINT ={str(round(this_LGINT, 2)).replace('-','−')}"]
                    ax.legend(handles, labels, loc='lower left', handlelength=0, handletextpad=0)
                    # ax.text(0.05, 0.85, f"$E_U =${this_EU_K:0.2f}", transform=ax.transAxes)
                    # ax.text(0.05, 0.75, f"$log_{{10}}(A_{{i,j}}) =${str(round(this_aij, 2)).replace('-','−')}", 
                    #         transform=ax.transAxes)
                    # if ret_tbl:
                    #     ax.text(0.05, 0.75, f"{tbl['Ju'][-1]}({tbl['Ku'][-1]}) – {tbl['Jl'][-1]}({tbl['Kl'][-1]})", 
                    #             transform=ax.transAxes)
                    # Show the local 2-sigma mark on each axis if requested (for upper-limit measurements)
                    if show_2_sigma:
                        # Get noise level for spectrum
                        noise_level = mad_std(data_sp_contsub.data)*u.K
                        ax.hlines(y=2*noise_level.value, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], 
                                  linewidth=1, color='r', linestyle='--')
                    # Record EUs and N_uppers if requested
                    if calc_N_uppers:
                        EUs.append(this_EU_K)
                        spec = pyspeckit.Spectrum(data=data_sp_contsub, 
                                                  xarr=molcube_slice.spectral_axis.to(u.GHz), unit=u.K)
                        spec.xarr.convert_to_unit(u.km/u.s, refX = freq, velocity_convention = 'radio')
                        spec_slice = spec.slice(v_cen - 4.5*u.km/u.s, v_cen + 4.5*u.km/u.s)
                        # Calculate integrated intensity (sum slice and multiply by channel width)
                        spec_slice_sum = spec_slice.data.sum()*u.K
                        channel_width = np.abs(spec_slice.xarr.cdelt())
                        mom0 = spec_slice_sum*channel_width
                        # Calculate upper state column density from integrated intensity, convert to logscale
                        N_upper = nupper_of_kkms(mom0, freq, 10**this_aij)
                        log_N_upper_g = np.log10(N_upper.value/this_deg)
                        log_N_upper_gs.append(log_N_upper_g)
                    if return_freqs:
                        mom0_freqs.append(freq)
                    ind+=1
    # Save the spectrum grid figures
    plt.tight_layout()
    name_for_file = name_for_plot.replace('$','').replace('^','').replace('{','').replace('}','').replace('_','')
    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/molecule_table/spectra_grids/grid_'+name_for_file+'_'+extr_type+'.pdf',
                facecolor = 'w', edgecolor = 'w', bbox_inches = 'tight')
    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/molecule_table/spectra_grids/grid_'+name_for_file+'_'+extr_type+'.png', 
                dpi = 250, bbox_inches = 'tight')
    # Export the upper state column density if requested
    if calc_N_uppers:
        return EUs, log_N_upper_gs
    # Export the frequencies for moment maps if requested
    if return_freqs:
        return mom0_freqs

def plot_mom0s(cubes, freqs, v_cen, fig_width, fig_height, nrows, ncols, name_for_plot=None, reg=None, 
               y=0.98, v_half=4.5):
    freqs_vals = [freq.value for freq in freqs]
    freqs = np.unique(np.array(freqs_vals))*u.MHz
    fig = plt.figure(figsize = (fig_width, fig_height))
    plt.suptitle(f"{name_for_plot}, $v_{{cen}} =$ {v_cen.value} km s$^{{-1}}$", fontsize = 16, y=y)
    ind=0 # Need this for subplots
    # Loop through cubes and extract lines
    for molcube in cubes:
        for freq in freqs:
            if molcube.spectral_axis.min() <= freq <= molcube.spectral_axis.max():
                ax = plt.subplot(nrows, ncols, ind+1)
                # Make moment map of line
                # molcube = molcube - molcube.median(axis=0, iterate_rays=True) # Continuum subtraction
                molcube_slice = molcube.with_spectral_unit(unit=u.km/u.s, velocity_convention='radio', 
                                                           rest_value=freq).spectral_slab(v_cen - v_half*u.km/u.s, v_cen + v_half*u.km/u.s)
                mom0 = molcube_slice.moment0()
                ax.imshow(mom0.value, origin='lower', cmap='inferno')
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
                if reg is not None:
                    ax.contour(reg, levels = [0, 1], linewidths=0.75, colors = ['c'])
                ax.text(0.05, 0.92, f"{freq.to(u.GHz):0.2f}", transform=ax.transAxes)
                ind+=1
    # Save the moment 0 map grid figures
    plt.tight_layout()
    name_for_file = name_for_plot.replace('$','').replace('^','').replace('{','').replace('}','').replace('_','')
    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/molecule_table/moment0_map_grids/mgrid_'+name_for_file+'_'+str(v_cen.value)+'.pdf', 
                facecolor = 'w', edgecolor = 'w', bbox_inches = 'tight')
    plt.savefig('/blue/adamginsburg/abulatek/brick/first_results/molecule_table/moment0_map_grids/mgrid_'+name_for_file+'_'+str(v_cen.value)+'.png', 
                dpi = 250, bbox_inches = 'tight')
    
def list_mol_tags(catalog):
    if catalog == 'JPL':
        from astroquery.jplspec import JPLSpec as QueryTool
    elif catalog == 'CDMS':
        from astroquery.linelists.cdms import CDMS as QueryTool
    else:
        raise ValueError("Invalid catalog specification")
    speciestab = QueryTool.get_species_table()
    speciestab.pprint(max_lines=-1)
    return speciestab

# def get_molecular_parameters_plus(catalog, fmin, fmax, mol_name_pyspeckit=None, mol_tag_pyspeckit=None):
#     # This is unnecessary because Adam pushed changes to get_molecular_parameters that fixed it; no need for hacking
#     molecule_name = mol_name_pyspeckit
#     molecule_tag = mol_tag_pyspeckit
#     if catalog == 'JPL':
#         from astroquery.jplspec import JPLSpec as QueryTool
#     elif catalog == 'CDMS':
#         from astroquery.linelists.cdms import CDMS as QueryTool
#     else:
#         raise ValueError("Invalid catalog specification")
#     speciestab = QueryTool.get_species_table()
#     if 'NAME' in speciestab.colnames:
#         molcol = 'NAME'
#     elif 'molecule' in speciestab.colnames:
#         molcol = 'molecule'
#     else:
#         raise ValueError(f"Did not find NAME or molecule in table columns: {speciestab.colnames}")
#     if molecule_tag is not None:
#         tagcol = 'tag' if 'tag' in speciestab.colnames else 'TAG'
#         match = speciestab[tagcol] == molecule_tag
#         molecule_name = speciestab[match][molcol][0]
#         if catalog == 'CDMS':
#             molsearchname = f'{molecule_tag:06d} {molecule_name}'
#         else:
#             molsearchname = f'{molecule_tag} {molecule_name}'
#         parse_names_locally = False
#         if molecule_name is not None:
#             print(f"WARNING: molecule_tag overrides molecule_name.  New molecule_name={molecule_name}.  Searchname = {molsearchname}")
#         else:
#             print(f"molecule_name={molecule_name} for tag molecule_tag={molecule_tag}.  Searchname = {molsearchname}")
#     else:
#         molsearchname = molecule_name
#         match = speciestab[molcol] == molecule_name
#         if match.sum() == 0:
#             # retry using partial string matching
#             match = np.core.defchararray.find(speciestab[molcol], molecule_name) != -1
#     if match.sum() != 1:
#         raise ValueError(f"Too many or too few matches ({match.sum()}) to {molecule_name}")
#     jpltable = speciestab[match]
#     jpltbl = QueryTool.query_lines(fmin, fmax, molecule=molsearchname,
#                                    parse_name_locally=parse_name_locally)
#     return jpltbl