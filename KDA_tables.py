import os
import numpy as np
from astropy.table import Table
from astropy import units as u


def get_kda_columns(t, keep_indices, colname='', dcolname='', factor=1.):
    col = t[colname][keep_indices]
    dcol = t[dcolname][keep_indices] * factor
    colmin = col - dcol
    colmax = col + dcol
    return col, dcol, colmin, colmax


def kda_info_table(t, glon='', glat='', vlsr='', dvlsr='', kda='',
                   a='', b='', pa='', factor_dvlsr=1., factor_a=1., factor_b=1.,
                   unit_dvlsr=u.km/u.s, unit_a=u.deg, unit_b=u.deg,
                   unit_pa=u.deg):
    p_far = t[kda].copy()
    keep = np.where(np.logical_or(p_far == 'N', p_far == 'F'))[0]

    p_far = np.where(p_far == 'N', -0.5, p_far)
    p_far = np.where(p_far == 'F', 0.5, p_far)
    p_far = p_far[keep].astype('float')

    factor_a = factor_a * (unit_a).to(u.deg)
    factor_b = factor_b * (unit_b).to(u.deg)
    dglon = dglat = a
    factor_dglon = factor_dglat = factor_a
    if pa:
        col_pa = np.deg2rad(t[pa][keep]) + np.pi / 2
    else:
        col_pa = np.zeros(len(keep))

    added_columns = [np.cos(col_pa), np.sin(col_pa),
                     (t[a][keep] * factor_a)**2, (t[b][keep] * factor_b)**2]
    added_names = ['cos_pa', 'sin_pa', 'aa', 'bb']

    factor_dvlsr = factor_dvlsr * (unit_dvlsr).to(u.km / u.s)

    col_glon, col_dglon, col_glon_min, col_glon_max = get_kda_columns(
        t, keep, colname=glon, dcolname=dglon, factor=factor_dglon)
    col_glat, col_dglat, col_glat_min, col_glat_max = get_kda_columns(
        t, keep, colname=glat, dcolname=dglat, factor=factor_dglat)
    col_vlsr, col_dvlsr, col_vlsr_min, col_vlsr_max = get_kda_columns(
        t, keep, colname=vlsr, dcolname=dvlsr, factor=factor_dvlsr)
    col_a = factor_a * t[a][keep]
    col_b = factor_a * t[b][keep]

    data = [col_glon, col_glat, col_vlsr, col_a, col_b, col_pa, col_dvlsr, p_far] + added_columns
    names = ['GLON', 'GLAT', 'VLSR', 'a', 'b', 'pa', 'd_VLSR', 'p_far'] + added_names

    t_sel = Table(data=data, names=names)

    for col in ['GLON', 'GLAT', 'a', 'b', 'pa']:
        t_sel[col].unit = 'deg'
        t_sel[col].format = "{0:.4f}"
    for col in ['VLSR', 'd_VLSR']:
        t_sel[col].unit = 'km/s'
        t_sel[col].format = "{0:.4f}"

    return t_sel
