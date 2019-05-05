#!/usr/bin/env python
"""
Plots of Fortran Venus atmosphere model
"""
import numpy as np
import xarray
from aeroplanets import venus

try:
    from matplotlib.pyplot import figure, show
except Exception:
    show = None


def plot_altprofile(atmo: xarray.Dataset):

    fg = figure(1, clear=True)
    ax1, ax2 = fg.subplots(1, 2, sharey=True)

    for v in atmo:
        if not v.startswith('T'):
            ax1.semilogx(atmo[v], atmo.alt_km, label=v)
    ax1.set_xlabel('density [cm$^{-3}$]')
    ax1.set_ylabel('altitude [km]')
    ax1.set_title(f'Venus VTS3 model: Density')
    ax1.set_xlim((1, 1e25))
    ax1.set_xscale('log')
    ax1.legend(loc='best')
    ax1.grid(True)

    for v in atmo:
        if v.startswith('T'):
            ax2.plot(atmo[v], atmo.alt_km, label=v)
    ax2.set_xlabel('temperature [K]')
    ax2.set_title(f'Temperature')
    ax2.legend(loc='best')
    ax2.grid(True)

    fg.suptitle(f'Venus VTS3 Atmosphere model: Hour {atmo.hour}'
                f' latitude: {atmo.lat} \nf107a, f107: {atmo.f107a}, {atmo.f107}')


if __name__ == '__main__':
    hour = 3.5
    alt_km = np.arange(1., 500., 1.)
    lat = 35.
    f107a = 100.
    f107 = 78.

    atmo = venus.venus_altitude(alt_km, lat, hour, f107a, f107)

    if atmo.alt_km.size > 1 and show is not None:
        plot_altprofile(atmo)
        show()
    else:
        print('Venus, local hour', atmo.hour, 'alt [km]:', atmo.alt_km)
        print('density: ', atmo.dens)
        print('temperatures:', atmo.temp)
