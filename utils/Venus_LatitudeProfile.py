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


def plot_latprofile(atmo: xarray.Dataset):

    fg = figure(1, clear=True)
    ax1, ax2 = fg.subplots(2, 1, sharex=True)

    for v in atmo:
        if not v.startswith('T'):
            ax1.semilogy(atmo.lat_deg, atmo[v], label=v)
    ax1.set_ylabel('density [cm$^{-3}$]')
    ax1.set_title(f'Density')
    ax1.legend(loc='best')
    ax1.grid(True)

    for v in atmo:
        if v.startswith('T'):
            ax2.plot(atmo.lat_deg, atmo[v], label=v)
    ax2.set_ylabel('temperature [K]')
    ax2.set_title(f'Temperature')
    ax2.legend(loc='best')
    ax2.grid(True)

    ax2.set_xlabel('Latitude [deg]')

    ax2.set_xlim((-90, 90))

    fg.suptitle(f'Venus VTS3 Atmosphere model: Altitude [km] {atmo.alt_km}'
                f' Hour {atmo.hour} \nf107a, f107: {atmo.f107a}, {atmo.f107}')


if __name__ == '__main__':
    hour = 3.5
    alt_km = 100.
    lat = np.arange(-90., 90. + 1., 1.0)
    f107a = 100.
    f107 = 78.

    atmo = venus.venus_latitude(alt_km, lat, hour, f107a, f107)

    if atmo.lat_deg.size > 1 and show is not None:
        plot_latprofile(atmo)
        show()
    else:
        print('Venus, local hour', atmo.hour, 'alt [km]:', atmo.alt_km)
        print('density: ', atmo.dens)
        print('temperatures:', atmo.temp)
