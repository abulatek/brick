# https://github.com/keflavich/Orion_ALMA_2016.1.00165.S/blob/master/analysis/lte_modeling_tools.py

import numpy as np
from astropy import constants
from astropy.modeling.models import custom_model
from astropy import units as u

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
        Nu (Q_rot / gu) exp(E_u/k Tex) = (3 k) / (8 pi^3 nu S mu^2 R_i)   (Q/g) exp(E_u/k Tex) integ(T_R/f dv)
    This term cancels:
        (Q_rot / gu) exp(E_u/k Tex)
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

def kkms_of_nupper(nupper, freq, Aul):
    """
    Convert the column density in the upper state of a line ``nupper`` to the
    integrated intensity in brightness units (K km / s).
    Inversion of nupper_of_kkms above.
    """

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    nupper = u.Quantity(nupper, u.cm**-2)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    Khz = (nupper / nline)

    kkms = (Khz / (freq/constants.c)).to(u.K * u.km/u.s)

    return kkms


def rovib_lte_model_generator(vibenergies, rotenergies):

    @custom_model
    def model(jstate, vstate, logcolumn=np.log(1e13), rottem=100, vibtem=2000):
        elower_vib = np.array([vibenergies[int(v)] for v in vstate])
        eupper_j = np.array([rotenergies[ju] for ju in jstate])
        #result = -1/rottem * (eupper - elower_vib) + column - 1/vibtem * eupper

        # these are the populations of states in the v=0 state at a given
        # J-level.
        result_v0 = -1/rottem * eupper_j + logcolumn

        # Then, for each of these, we determine the levels in the vibrationally
        # excited state by adding e^-hnu/kt, where t is a different temperature
        # (t_v) and nu is now just the nu provided by the vibrations
        result = result_v0 - 1/vibtem * elower_vib

        return result

    return model()

def simple_lte_model_generator():

    @custom_model
    def model(eupper, logcolumn=np.log(1e13), tem=100):
        """
        Calculate the quantity N_u/g_u as a function of E_u in Kelvin
        The 'logcolumn' quantity is N_tot / Q_tot
        Temperature is the excitation temperature
        """

        result = -1/tem * eupper + logcolumn

        return result

    return model()

def get_molecular_parameters(molecule_name, tex=50, fmin=1*u.GHz, fmax=1*u.THz,
                             catalog='JPL',**kwargs):
    """
    Get the molecular parameters for a molecule from the JPL or CDMS catalog
    (this version should, in principle, be entirely self-consistent)
    Parameters
    ----------
    molecule_name : string
        The string name of the molecule (normal name, like CH3OH or CH3CH2OH,
        but it has to match the JPL catalog spec)
    tex : float
        Optional excitation temperature (basically checks if the partition
        function calculator works)
    catalog : 'JPL' or 'CDMS'
        Which catalog to pull from
    fmin : quantity with frequency units
    fmax : quantity with frequency units
        The minimum and maximum frequency to search over
    Examples
    --------
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH2CHCN',
    ...                                                          fmin=220*u.GHz,
    ...                                                          fmax=222*u.GHz,
                                                                )
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
    ...                                                          fmin=90*u.GHz,
    ...                                                          fmax=100*u.GHz)
    """
    if catalog == 'JPL':
        from astroquery.jplspec import JPLSpec as QueryTool
    elif catalog == 'CDMS':
        from astroquery.cdms import CDMS as QueryTool
    else:
        raise ValueError("Invalid catalog specification")

    speciestab = QueryTool.get_species_table()
    jpltable = speciestab[speciestab['NAME'] == molecule_name]
    if len(jpltable) != 1:
        raise ValueError(f"Too many or too few matches to {molecule_name}")

    jpltbl = QueryTool.query_lines(fmin, fmax, molecule=molecule_name,
                                   parse_name_locally=True)

    def partfunc(tem):
        """
        interpolate the partition function
        WARNING: this can be very wrong
        """
        tem = u.Quantity(tem, u.K).value
        tems = np.array(jpltable.meta['Temperature (K)'])
        keys = [k for k in jpltable.keys() if 'q' in k.lower()]
        logQs = jpltable[keys]
        logQs = np.array(list(logQs[0]))
        inds = np.argsort(tems)
        #logQ = np.interp(tem, tems[inds], logQs[inds])
        # linear interpolation is appropriate; Q is linear with T... for some cases...
        # it's a safer interpolation, anyway.
        # to get a better solution, you can fit a functional form as shown in the
        # JPLSpec docs, but that is... left as an exercise.
        # (we can test the degree of deviation there)
        linQ = np.interp(tem, tems[inds], 10**logQs[inds])
        return linQ

    freqs = jpltbl['FREQ'].quantity
    freq_MHz = freqs.to(u.MHz).value
    deg = np.array(jpltbl['GUP'])
    EL = jpltbl['ELO'].quantity.to(u.erg, u.spectral())
    dE = freqs.to(u.erg, u.spectral())
    EU = EL + dE

    # need elower, eupper in inverse centimeter units
    elower_icm = jpltbl['ELO'].quantity.to(u.cm**-1).value
    eupper_icm = elower_icm + (freqs.to(u.cm**-1, u.spectral()).value)

    # from Brett McGuire https://github.com/bmcguir2/simulate_lte/blob/1f3f7c666946bc88c8d83c53389556a4c75c2bbd/simulate_lte.py#L2580-L2587

    # LGINT: Base 10 logarithm of the integrated intensity in units of nm2 ·MHz at 300 K.
    # (See Section 3 for conversions to other units.)
    # see also https://cdms.astro.uni-koeln.de/classic/predictions/description.html#description
    CT = 300
    logint = np.array(jpltbl['LGINT']) # this should just be a number
    #from CDMS website
    sijmu = (np.exp(np.float64(-(elower_icm/0.695)/CT)) - np.exp(np.float64(-(eupper_icm/0.695)/CT)))**(-1) * ((10**logint)/freq_MHz) * (24025.120666) * partfunc(CT)

    #aij formula from CDMS.  Verfied it matched spalatalogue's values
    aij = 1.16395 * 10**(-20) * freq_MHz**3 * sijmu / deg

    # we want logA for consistency with use in generate_model below
    aij = np.log10(aij)
    EU = EU.to(u.erg).value

    ok = np.isfinite(aij) & np.isfinite(EU) & np.isfinite(deg) & np.isfinite(freqs)

    return freqs[ok], aij[ok], deg[ok], EU[ok], partfunc # What is partfunc supposed to be? I want it for each line, not the function itself



if __name__ == '__main__':
    # round-trip test
    kkms = 100*u.K*u.km/u.s
    freq = 100*u.GHz
    Aul = 1*u.s**-1
    degeneracies = 1
    nupper = nupper_of_kkms(kkms, freq, Aul, degeneracies)
    kkms2 = kkms_of_nupper(nupper, freq, Aul, degeneracies)
    np.testing.assert_almost_equal(kkms2.value, kkms.value)