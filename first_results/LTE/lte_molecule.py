"""
LTE Molecule Modeling Tool
==========================

Uses astroquery & vamdclib to obtain molecular parameters.
http://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html

Equations are based on Mangum & Shirley 2015 (2015PASP..127..266M)

Module API
^^^^^^^^^^
"""
from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import constants
from .model import SpectralModel

kb_cgs = constants.k_B.cgs.value
h_cgs = constants.h.cgs.value
eightpicubed = 8 * np.pi**3
threehc = 3 * constants.h.cgs * constants.c.cgs
hoverk_cgs = (h_cgs/kb_cgs)
c_cgs = constants.c.cgs.value

def line_tau(tex, total_column, partition_function, degeneracy, frequency,
             energy_upper, einstein_A=None):
    """
    Given the excitation temperature of the state, total column density of the
    molecule, the partition function, the degeneracy of the state, the
    frequency of the state, and the upper-state energy level, return the optical
    depth of that transition.

    This is a helper function for the LTE molecule calculations.  It implements
    the equations

    .. math::

        \\int \\tau_\\nu d\\nu = \\frac{c^2}{8 \pi \\nu^2} A_{ij} N_u
                    \\left[ \\exp\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    or

    .. math::

        \\tau_\\nu = \\frac{c^2}{8 \pi \\nu^2} A_{ij} N_u \\phi_\\nu
                    \\left[ \\exp\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    where

    .. math::
        N_{u} = N_{tot} \\frac{g_u}{Q} \\exp\left(\\frac{-E_u}{k_B T_{ex}} \\right)

    based on Equations 11 and 29 of Mangum & Shirley 2015 (2015PASP..127..266M)

    The return value is therefore

    .. math::

        \\tau_\\nu / \\phi_\\nu


    The line profile function is, for a Gaussian, given by eqn A1:

    .. math::

        \\phi_\\nu = \\frac{1}{\\sqrt{2 \\pi} \\sigma}
        \\exp \\left[ -\\frac{(\\nu-\\nu_0)^2}{2 \\sigma^2} \\right]

    """
    # don't use dipole moment, because there are extra hidden dependencies

    assert frequency.unit.is_equivalent(u.Hz)
    assert energy_upper.unit.is_equivalent(u.erg)
    assert total_column.unit.is_equivalent(u.cm**-2)
    assert tex.unit.is_equivalent(u.K)

    N_upper = (total_column * degeneracy / partition_function *
               np.exp(-energy_upper / (constants.k_B * tex)))

    # equation 29 in Mangum 2015
    #taudnu = (eightpicubed * frequency * dipole_moment**2 / threehc *
    #             (np.exp(frequency*h_cgs/(kb_cgs*tex))-1) * N_upper)
    assert einstein_A.unit.is_equivalent(u.Hz)
    taudnu = ((constants.c**2/(8*np.pi*frequency**2) * einstein_A * N_upper)*
              (np.exp(frequency*constants.h/(constants.k_B*tex))-1))

    return taudnu.decompose()

# Deprecated version of the above
# def line_tau_nonquantum(tex, total_column, partition_function, degeneracy,
#                         frequency, energy_upper, SijMu2=None, molwt=None):
#
#     assert frequency.unit.is_equivalent(u.Hz)
#     assert energy_upper.unit.is_equivalent(u.erg)
#     assert total_column.unit.is_equivalent(u.cm**-2)
#     assert tex.unit.is_equivalent(u.K)
#
#     energy_lower = energy_upper - frequency*constants.h
#     #N_lower = (total_column * degeneracy / partition_function *
#     #           np.exp(-energy_lower / (constants.k_B * tex)))
#
#     # http://www.cv.nrao.edu/php/splat/OSU_Splat.html
#     assert SijMu2.unit.is_equivalent(u.debye**2)
#     amu = u.Da
#     assert molwt.unit.is_equivalent(amu)
#     C1 = 54.5953 * u.nm**2 * u.K**0.5 / amu**0.5 / u.debye**2
#     C2 = 4.799237e-5 * u.K / u.MHz
#     C3 = (1.43877506 * u.K / ((1*u.cm).to(u.Hz, u.spectral()) * constants.h)).to(u.K/u.erg)
#     #C3 = (constants.h / constants.k_B).to(u.K/u.erg)
#     tau = total_column/partition_function * C1 * (molwt/tex)**0.5 * (1-np.exp(-C2*frequency/tex)) * SijMu2 * np.exp(-C3*energy_lower/tex)
#
#     return tau.decompose()

def line_tau_cgs(tex, total_column, partition_function, degeneracy, frequency,
                 energy_upper, einstein_A):
    """
    Given the excitation temperature of the state, total column density of the
    molecule, the partition function, the degeneracy of the state, the
    frequency of the state, and the upper-state energy level, return the optical
    depth of that transition.

    Unlike :func:`line_tau`, this function requires inputs in CGS units.

    This is a helper function for the LTE molecule calculations.  It implements
    the equations

    .. math::

        \\int \\tau_\\nu d\\nu = \\frac{c^2}{8 \pi \\nu^2} A_{ij} N_u
                    \\left[ \\exp\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    or

    .. math::

        \\tau_\\nu = \\frac{c^2}{8 \pi \\nu^2} A_{ij} N_u \\phi_\\nu
                    \\left[ \\exp\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    where

    .. math::
        N_{u} = N_{tot} \\frac{g_u}{Q} \\exp\left(\\frac{-E_u}{k_B T_{ex}} \\right)

    based on Equations 11 and 29 of Mangum & Shirley 2015 (2015PASP..127..266M)

    The return value is therefore

    .. math::

        \\tau_\\nu / \\phi_\\nu


    The line profile function is, for a Gaussian, given by eqn A1:

    .. math::

        \\phi_\\nu = \\frac{1}{\\sqrt{2 \\pi} \\sigma}
        \\exp \\left[ -\\frac{(\\nu-\\nu_0)^2}{2 \\sigma^2} \\right]

    """

    N_upper = (total_column * degeneracy / partition_function *
               np.exp(-energy_upper / (kb_cgs * tex)))

    # equation 29 in Mangum 2015
    #taudnu = (eightpicubed * frequency * dipole_moment**2 / threehc *
    #             (np.exp(frequency*h_cgs/(kb_cgs*tex))-1) * N_upper)

    # substitute eqn 11:
    # Aij = 64 pi^4 nu^3 / (3 h c^3) |mu ul|^2
    # becomes
    # dipole_moment^2 = 3 h c^3 Aij / ( 64 pi^4 nu^3 )

    taudnu = ((c_cgs**2/(8*np.pi*frequency**2) * einstein_A * N_upper)*
              (np.exp(frequency*h_cgs/(kb_cgs*tex))-1))

    return taudnu

def Jnu(nu, T):
    """RJ equivalent temperature (MS15 eqn 24)"""
    return constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*T))-1)

def Jnu_cgs(nu, T):
    """RJ equivalent temperature (MS15 eqn 24)
    (use cgs constants for speed)
    """
    return hoverk_cgs*nu / (np.exp(hoverk_cgs*nu/T)-1)

def line_brightness(tex, dnu, frequency, tbg=2.73*u.K, *args, **kwargs):
    tau = line_tau(tex=tex, frequency=frequency, *args, **kwargs) / dnu
    tau = tau.decompose()
    assert tau.unit == u.dimensionless_unscaled
    return (Jnu(frequency, tex)-Jnu(frequency, tbg)).decompose() * (1 - np.exp(-tau))

def line_brightness_cgs(tex, dnu, frequency, tbg=2.73, *args, **kwargs):
    tau = line_tau(tex=tex, frequency=frequency, *args, **kwargs) / dnu
    return (Jnu(frequency, tex)-Jnu(frequency, tbg)) * (1 - np.exp(-tau))

# requires vamdc branch of astroquery
def get_molecular_parameters(molecule_name,
                             molecule_name_jpl=None,
                             tex=50, fmin=1*u.GHz, fmax=1*u.THz,
                             line_lists=['SLAIM'],
                             export_limit=1e5,
                             chem_re_flags=0, **kwargs):
    """
    Get the molecular parameters for a molecule from the CDMS database using
    vamdclib

    Parameters
    ----------
    molecule_name : string
        The string name of the molecule (normal name, like CH3OH or CH3CH2OH)
    molecule_name_jpl : string or None
        Default to molecule_name, but if jplspec doesn't have your molecule
        by the right name, specify it here
    tex : float
        Optional excitation temperature (basically checks if the partition
        function calculator works)
    fmin : quantity with frequency units
    fmax : quantity with frequency units
        The minimum and maximum frequency to search over
    line_lists : list
        A list of Splatalogue line list catalogs to search.  Valid options
        include SLAIM, CDMS, JPL.  Only a single catalog should be used to
        avoid repetition of transitions and species
    chem_re_flags : int
        An integer flag to be passed to splatalogue's chemical name matching
        tool

    Examples
    --------
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters(molecule_name='CH2CHCN',
    ...                                                          fmin=220*u.GHz,
    ...                                                          fmax=222*u.GHz,
                                                                )
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
    ...                                                          fmin=90*u.GHz,
    ...                                                          fmax=100*u.GHz)
    """
    from astroquery.splatalogue import Splatalogue
    from astroquery.jplspec import JPLSpec

    # do this here, before trying to compute the partition function, because
    # this query could fail
    tbl = Splatalogue.query_lines(fmin, fmax, chemical_name=molecule_name,
                                  line_lists=line_lists,
                                  show_upper_degeneracy=True,
                                  export_limit=export_limit,
                                  **kwargs)

    if molecule_name_jpl is None:
        molecule_name_jpl = molecule_name

    jpltable = JPLSpec.get_species_table()[JPLSpec.get_species_table()['NAME'] == molecule_name_jpl]
    if len(jpltable) != 1:
        raise ValueError(f"Too many or too few matches to {molecule_name_jpl}")

    def partfunc(tem):
        """
        interpolate the partition function
        WARNING: this can be very wrong
        """
        tem = u.Quantity(tem, u.K).value
        tems = np.array(jpltable.meta['Temperature (K)'])
        logQs = jpltable['QLOG1 QLOG2 QLOG3 QLOG4 QLOG5 QLOG6 QLOG7'.split()]
        logQs = np.array(list(logQs[0]))
        inds = np.argsort(tems)
        logQ = np.interp(tem, tems[inds], logQs[inds])
        return 10**logQ


    freqs = (np.array(tbl['Freq-GHz'])*u.GHz if 'Freq-GHz' in tbl.colnames else
             np.array(tbl['Freq-GHz(rest frame,redshifted)'])*u.GHz)
    aij = tbl['Log<sub>10</sub> (A<sub>ij</sub>)']
    deg = tbl['Upper State Degeneracy']
    EU = (np.array(tbl['E_U (K)'])*u.K*constants.k_B).to(u.erg).value

    return freqs, aij, deg, EU, partfunc


def get_molecular_parameters_JPL(molecule_name_jpl, tex=50, fmin=1*u.GHz,
                                 fmax=1*u.THz, **kwargs):
    """
    Get the molecular parameters for a molecule from the JPL catalog

    (this version should, in principle, be entirely self-consistent)

    Parameters
    ----------
    molecule_name : string
        The string name of the molecule (normal name, like CH3OH or CH3CH2OH,
        but it has to match the JPL catalog spec)
    tex : float
        Optional excitation temperature (basically checks if the partition
        function calculator works)
    fmin : quantity with frequency units
    fmax : quantity with frequency units
        The minimum and maximum frequency to search over

    Examples
    --------
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters(molecule_name='CH2CHCN',
    ...                                                          fmin=220*u.GHz,
    ...                                                          fmax=222*u.GHz,
                                                                )
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
    ...                                                          fmin=90*u.GHz,
    ...                                                          fmax=100*u.GHz)
    """
    from astroquery.jplspec import JPLSpec

    jpltable = JPLSpec.get_species_table()[JPLSpec.get_species_table()['NAME'] == molecule_name_jpl]
    if len(jpltable) != 1:
        raise ValueError(f"Too many or too few matches to {molecule_name_jpl}")

    jpltbl = JPLSpec.query_lines(fmin, fmax, molecule=molecule_name_jpl,
                                 parse_name_locally=True)

    def partfunc(tem):
        """
        interpolate the partition function
        WARNING: this can be very wrong
        """
        tem = u.Quantity(tem, u.K).value
        tems = np.array(jpltable.meta['Temperature (K)'])
        logQs = jpltable['QLOG1 QLOG2 QLOG3 QLOG4 QLOG5 QLOG6 QLOG7'.split()]
        logQs = np.array(list(logQs[0]))
        inds = np.argsort(tems)
        logQ = np.interp(tem, tems[inds], logQs[inds])
        return 10**logQ

    freqs = jpltbl['FREQ'].quantity
    freq_MHz = freqs.to(u.MHz).value
    deg = jpltbl['GUP'].value
    EL = jpltbl['ELO'].quantity.to(u.erg, u.spectral())
    dE = freqs.to(u.erg, u.spectral())
    EU = EL + dE

    # need elower, eupper in inverse centimeter units
    elower_icm = jpltbl['ELO'].quantity.to(u.cm**-1).value
    eupper_icm = elower_icm + (freqs.to(u.cm**-1, u.spectral()).value)

    # from Brett McGuire https://github.com/bmcguir2/simulate_lte/blob/1f3f7c666946bc88c8d83c53389556a4c75c2bbd/simulate_lte.py#L2580-L2587

    # LGINT: Base 10 logarithm of the integrated intensity in units of nm2 ·MHz at 300 K. (See Section 3 for conversions to other units.)
    CT = 300
    logint = jpltbl['LGINT'].value
    #from CDMS website
    sijmu = (np.exp(np.float64(-(elower_icm/0.695)/CT)) - np.exp(np.float64(-(eupper_icm/0.695)/CT)))**(-1) * ((10**logint)/freq_MHz) * (24025.120666) * partfunc(CT)

    #aij formula from CDMS.  Verfied it matched spalatalogue's values
    aij = 1.16395 * 10**(-20) * freq_MHz**3 * sijmu / deg

    # we want logA for consistency with use in generate_model below
    aij = np.log10(aij)
    EU = EU.to(u.erg).value

    return freqs, aij, deg, EU, partfunc



def generate_model(xarr, vcen, width, tex, column,
                   freqs, aij, deg, EU, partfunc,
                   background=None, tbg=2.73,
                  ):
    """
    Model Generator
    """

    if hasattr(tex,'unit'):
        tex = tex.value
    if hasattr(tbg,'unit'):
        tbg = tbg.value
    if hasattr(column, 'unit'):
        column = column.value
    if column < 25:
        column = 10**column
    if hasattr(vcen, 'unit'):
        vcen = vcen.value
    if hasattr(width, 'unit'):
        width = width.value

    ckms = constants.c.to(u.km/u.s).value

    # assume equal-width channels
    #kwargs = dict(rest=ref_freq)
    #equiv = u.doppler_radio(**kwargs)

    # channelwidth array, with last element approximated
    channelwidth = np.empty_like(xarr.value)
    channelwidth[:-1] = np.abs(np.diff(xarr.to(u.Hz))).value
    channelwidth[-1] = channelwidth[-2]

    #velo = xarr.to(u.km/u.s, equiv).value
    freq = xarr.to(u.Hz).value # same unit as nu below
    model = np.zeros_like(xarr).value

    # splatalogue can report bad frequencies as zero
    OK = freqs.value != 0

    freqs_ = freqs.to(u.Hz).value

    Q = partfunc(tex)

    for logA, gg, restfreq, eu in zip(aij[OK], deg[OK], freqs_[OK], EU[OK]):
        tau_over_phi = line_tau_cgs(tex=tex, total_column=column,
                                    partition_function=Q, degeneracy=gg,
                                    frequency=restfreq, energy_upper=eu,
                                    einstein_A=10**logA)
        width_dnu = width / ckms * restfreq

        phi_nu = (
            ((2*np.pi)**0.5 * width_dnu)**-1 *
            np.exp(-(freq-(1-vcen/ckms)*restfreq)**2/(2*width_dnu**2)))

        tau_profile = (tau_over_phi * phi_nu)

        jnu = (Jnu_cgs(restfreq, tex)-Jnu_cgs(restfreq, tbg))

        model = model + jnu*(1-np.exp(-tau_profile))

    if background is not None:
        return background-model
    return model

"""
Example case to produce a model:

freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH')
def modfunc(xarr, vcen, width, tex, column):
    return generate_model(xarr, vcen, width, tex, column, freqs=freqs, aij=aij,
                          deg=deg, EU=EU, partfunc=partfunc)

fitter = generate_fitter(modfunc, name="CH3OH")
"""

def generate_fitter(model_func, name):
    """
    Generator for hnco fitter class
    """

    myclass = SpectralModel(model_func, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
    myclass.__name__ = name

    return myclass


def nupper_of_kkms(kkms, freq, Aul, replace_bad=None):
    """
    Mangum & Shirley 2015 eqn 82 gives, for the optically thin, Rayleigh-Jeans,
    negligible background approximation:

        Ntot = (3 k) / (8 pi^3 nu S mu^2 R_i)   (Q/g) exp(E_u/k Tex) integ(T_R/f dv)

    Eqn 31:

        Ntot/Nu = Q_rot / gu exp(E_u/k Tex)

        -> Ntot = Nu Q_rot / gu exp(E_u/k Tex)
        -> Nu = N_tot g / Qrot exp(-E_u / k Tex)

    To get Nu of an observed line, then:

        Nu Q_rot / gu exp(E_u/k Tex) = (3 k) / (8 pi^3 nu S mu^2 R_i)   (Q/g) exp(E_u/k Tex) integ(T_R/f dv)

    This term cancels:
        Q_rot / gu exp(E_u/k Tex)

    Leaving:

        Nu = (3 k) / (8 pi^3 nu S mu^2 R_i)   integ(T_R/f dv)

    integ(T_R/f dv) is the optically thin integrated intensity in K km/s
    dnu/nu = dv/c [doppler eqn], so to get integ(T_R dnu), sub in dv = c/nu dnu

        Nu = (3 k c) / (8 pi^3 nu^2  S mu^2 R_i)   integ(T_R/f dnu)


    We then need to deal with the S mu^2 R_i term.  We assume R_i = 1, since we
    are not measuring any hyperfine transitions (R_i is the hyperfine
    degeneracy; eqn 75)
    Equation 11:

        A_ul = 64 pi^4 nu^3 / (3 h c^3) |mu_ul|^2

    Equation 62:

        |mu_ul|^2 = S mu^2

        -> S mu^2 = (3 h c^3 Aul) / (64 pi^4 nu^3)

    Plugging that in gives

        Nu = (3 k c) / (8 pi^3 nu^2  ((3 h c^3 Aul) / (64 pi^4 nu^3)))   integ(T_R/f dnu)
           = (3 k c 64 pi^4 nu^3) / (8 pi^3 nu^2 3 h c^3 Aul)            integ(T_R/f dnu)
           = (8 pi nu k / (Aul c^2 h)) integ(T_R/f dnu)

    which is the equation implemented below.  We could also have left this in
    dv units by substituting du = nu/c dv:

           = (8 pi nu^2 k / (Aul c^3 h)) integ(T_R/f dv)

    """

    if replace_bad:
        neg = kkms <= 0
        kkms[neg] = replace_bad

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)

    return (nline * Khz).to(u.cm**-2)

def ntot_of_nupper(nupper, eupper, tex, Q_rot, degeneracy=1):
    """ Given an N_upper, E_upper, tex, Q_rot, and degeneracy for a single state, give N_tot

    Mangum & Shirley 2015 eqn 31
    Ntot/Nu = Q_rot / gu exp(E_u/k Tex)

    Example:
        >>> tex = 50*u.K
        >>> kkms = 100*u.K*u.km/u.s
        >>> from pyspeckit.spectrum.models import lte_molecule
        >>> freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(molecule_name='HNCO v=0', molecule_name_jpl='HNCO', fmin=87*u.GHz, fmax=88*u.GHz)
        >>> nupper = lte_molecule.nupper_of_kkms(kkms, freqs, 10**aij)
        >>> ntot = lte_molecule.ntot_of_nupper(nupper, EU*u.erg, tex, Q_rot=partfunc(tex), degeneracy=deg)
    """

    Ntot = nupper * (Q_rot/degeneracy) * np.exp(eupper / (constants.k_B*tex))

    return Ntot

# url = 'http://cdms.ph1.uni-koeln.de/cdms/tap/'
# rslt = requests.post(url+"/sync", data={'REQUEST':"doQuery", 'LANG': 'VSS2', 'FORMAT':'XSAMS', 'QUERY':"SELECT SPECIES WHERE MoleculeStoichiometricFormula='CH2O'"})

if __name__ == "__main__":
    # example

    # 303
    J = 3
    gI = 0.25
    gJ = 2*J+1
    gK = 1

    # 321 has same parameters for g

    ph2co = {'tex':18.75*u.K,
             'total_column': 1e12*u.cm**-2,
             'partition_function': 44.6812, # splatalogue's 18.75
             'degeneracy': gI*gJ*gK,
             #'dipole_moment': 2.331e-18*u.esu*u.cm, #2.331*u.debye,
            }

    ph2co_303 = {
             'frequency': 218.22219*u.GHz,
             'energy_upper': kb_cgs*20.95582*u.K,
             'einstein_A': 10**-3.55007/u.s,
    }
    ph2co_303.update(ph2co)
    ph2co_303['dnu'] = (1*u.km/u.s/constants.c * ph2co_303['frequency'])

    ph2co_321 = {
             'frequency': 218.76007*u.GHz,
             'energy_upper': kb_cgs*68.11081*u.K,
             'einstein_A': 10**-3.80235/u.s,
    }
    ph2co_321.update(ph2co)
    ph2co_321['dnu'] = (1*u.km/u.s/constants.c * ph2co_321['frequency'])

    ph2co_322 = {
             'frequency': 218.47563*u.GHz,
             'energy_upper': kb_cgs*68.0937*u.K,
             'einstein_A': 10**-3.80373/u.s,
    }
    ph2co_322.update(ph2co)
    ph2co_322['dnu'] = (1*u.km/u.s/constants.c * ph2co_322['frequency'])

    print(("tau303 = {0}".format(line_tau(**ph2co_303))))
    print(("tau321 = {0}".format(line_tau(**ph2co_321))))
    print(("tau322 = {0}".format(line_tau(**ph2co_322))))
    print(("r303/r321 = {0}".format(line_brightness(**ph2co_321)/line_brightness(**ph2co_303))))
    print(("r303/r322 = {0}".format(line_brightness(**ph2co_322)/line_brightness(**ph2co_303))))

    # CDMS Q
    import requests
    import bs4
    url = 'http://cdms.ph1.uni-koeln.de/cdms/tap/'
    rslt = requests.post(url+"/sync", data={'REQUEST':"doQuery", 'LANG': 'VSS2', 'FORMAT':'XSAMS', 'QUERY':"SELECT SPECIES WHERE MoleculeStoichiometricFormula='CH2O'"})
    bb = bs4.BeautifulSoup(rslt.content, 'html5lib')
    h = [x for x in bb.findAll('molecule') if x.ordinarystructuralformula.value.text=='H2CO'][0]
    tem_, Q_ = h.partitionfunction.findAll('datalist')
    tem = [float(x) for x in tem_.text.split()]
    Q = [float(x) for x in Q_.text.split()]

    del ph2co_303['tex']
    del ph2co_303['partition_function']
    T_303 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_303).value for tex,pf in
                      zip(tem,Q)])

    del ph2co_321['tex']
    del ph2co_321['partition_function']
    T_321 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_321).value for tex,pf in
                      zip(tem,Q)])

    del ph2co_322['tex']
    del ph2co_322['partition_function']
    T_322 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_322).value for tex,pf in
                      zip(tem,Q)])

    import pylab as pl

    pl.clf()
    pl.subplot(2,1,1)
    pl.plot(tem, T_321, label='$3_{2,1}-2_{2,0}$')
    pl.plot(tem, T_322, label='$3_{2,2}-2_{2,1}$')
    pl.plot(tem, T_303, label='$3_{0,3}-2_{0,2}$')
    pl.xlim(0,200)
    pl.subplot(2,1,2)
    pl.plot(tem, T_321/T_303, label='321/303')
    pl.plot(tem, T_322/T_303, label='322/303')
    pl.xlim(0,200)

    pl.draw(); pl.show()
