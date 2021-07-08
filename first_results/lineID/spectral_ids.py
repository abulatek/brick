from astroquery.splatalogue import Splatalogue, utils
from astropy import stats
from astropy import log
from astropy import constants
from astropy import units as u
from astropy import table
from astropy.table import Column
from astropy.io import fits
from scipy.ndimage import label
from radio_beam import Beams
import numpy as np

import pyspeckit
import glob

def jtok_factors(beams, freqs):
    factors = []
    for bm,frq in zip(beams, freqs):

        # create a beam equivalency for brightness temperature
        bmequiv = bm.jtok_equiv(frq)
        factor = (u.Jy).to(u.K, equivalencies=bmequiv)
        factors.append(factor)
    factor = np.array(factors)
    return factor

def line_ids(fn):
    results = []
    Ulines = []

    sp = pyspeckit.Spectrum(fn)

    # this is a bit hackier than I like
    # we'll do all our measurements in Kelvin!
    beams = Beams.from_fits_bintable(fits.open(fn)[1])
    factors = jtok_factors(beams, sp.xarr.to(u.GHz))
    sp.data = sp.data*factors
    sp.unit = u.K

    # want km/s - reference will be ~middle of SPW
    sp.xarr.convert_to_unit(u.km/u.s)

    med = np.nanmedian(sp.data)

    mad = stats.mad_std(sp.data - med)
    detections = (sp.data-med) > 5*mad

    labels, ct = label(detections)

    for labelid in range(1,ct+1):
        ssp = sp[labels == labelid]
        try:
            ssp.xarr.convert_to_unit(u.GHz)
            ssp.specfit()
            ssp.specfit.parinfo
            frq = ssp.specfit.parinfo['SHIFT0'].value * ssp.xarr.unit
        except Exception as ex:
            print(ex)
            frq = ssp.xarr.to(u.GHz).mean()
        sq = Splatalogue.query_lines(frq*(1+0/3e5), frq*(1+75/3e5), # 30/3e5 original lower bound
                                     only_astronomically_observed=True)
        if len(sq) > 0:
            tbl = utils.minimize_table(sq)
            try:
                total_intensity = ssp.data.sum() * np.abs(ssp.xarr.to(u.km/u.s).cdelt())
            except ValueError:
                total_intensity = ssp.data.sum() * np.abs(sp.xarr.to(u.km/u.s).cdelt())
            peak_intensity = ssp.data.max()
            tbl.add_column(Column(data=total_intensity, name='TotalIntensity'))
            tbl.add_column(Column(data=peak_intensity, name='PeakIntensity'))
            tbl.add_column(Column(data=mad, name='RMS'))
            tbl.add_column(Column(data=u.Quantity((-(frq.to(u.GHz).value -
                                                     tbl['Freq']) / tbl['Freq']
                                                   * constants.c), u.km/u.s),
                                  name='Velocity'))
#             print(tbl.pprint(max_width=200))
            results.append(tbl)
        else:
            log.warning(f"Frequency {frq.to(u.GHz)} had no hits")
            Ulines.append(frq)

    try:
        match_table = table.vstack(results)
    except ValueError:
        pass
    else:
#         match_table.remove_column('QNs')
        match_table = table.unique(match_table, keys='Species')
    match_table.sort('Freq')
    print(match_table.pprint(max_width=200))
    print(match_table['Species','Freq','QNs','Velocity'])
#     match_table.write(f"line_fit_table_{suffix}.ipac", format='ascii.ipac', overwrite=True)
    return match_table
