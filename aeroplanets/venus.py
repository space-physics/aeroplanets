import numpy as np
import numpy.f2py
import xarray
from pathlib import Path

try:
    import venus_vts3
except ImportError:
    from aeroplanets.builder import build
    build(Path(__file__).resolve().parents[1] / 'src/planet/Venus//pvatmos.f', 'venus_vts3')


def venus_altitude(alt_km: float, lat_deg: float, hour: float,
                   f107a: float, f107: float) -> xarray.Dataset:
    """
    Altitude profile of Venus neutral atmosphere
    """

    alt_km = np.atleast_1d(alt_km)
    dens = np.empty((alt_km.size, 7))
    temp = np.empty((alt_km.size, 2))
    for i in range(alt_km.size):
        dens[i, :], temp[i, :] = venus_vts3.vts3(alt_km[i], 0, lat_deg, hour, f107a, f107, 48)

    var = 'alt_km'
    atmo = xarray.Dataset(coords={var: alt_km})
    atmo['CO2'] = (var, dens[:, 1])
    atmo['O'] = (var, dens[:, 2])
    atmo['CO'] = (var, dens[:, 3])
    atmo['He'] = (var, dens[:, 4])
    atmo['N'] = (var, dens[:, 5])
    atmo['N2'] = (var, dens[:, 6])
    atmo['Tn'] = (var, temp[:, 1])

    atmo.attrs['lat'] = lat_deg
    atmo.attrs['f107a'] = f107a
    atmo.attrs['f107'] = f107
    atmo.attrs['hour'] = hour
    atmo.attrs['Texo'] = temp[0, 0]

    return atmo


def venus_time(alt_km: float, lat_deg: float, hour: float,
               f107a: float, f107: float) -> xarray.Dataset:
    """
    Time profile of Venus neutral atmosphere
    """

    hour = np.atleast_1d(hour)
    dens = np.empty((hour.size, 7))
    temp = np.empty((hour.size, 2))
    for i in range(hour.size):
        dens[i, :], temp[i, :] = venus_vts3.vts3(alt_km, 0, lat_deg, hour[i], f107a, f107, 48)

    var = 'hour'
    atmo = xarray.Dataset(coords={var: hour})
    atmo['CO2'] = (var, dens[:, 1])
    atmo['O'] = (var, dens[:, 2])
    atmo['CO'] = (var, dens[:, 3])
    atmo['He'] = (var, dens[:, 4])
    atmo['N'] = (var, dens[:, 5])
    atmo['N2'] = (var, dens[:, 6])

    atmo['Texo'] = (var, temp[:, 0])
    atmo['Tn'] = (var, temp[:, 1])

    atmo.attrs['lat'] = lat_deg
    atmo.attrs['f107a'] = f107a
    atmo.attrs['f107'] = f107
    atmo.attrs['alt_km'] = alt_km

    return atmo


def venus_latitude(alt_km: float, lat_deg: float, hour: float,
                   f107a: float, f107: float) -> xarray.Dataset:
    """
    Latitude profile of Venus neutral atmosphere
    """

    lat_deg = np.atleast_1d(lat_deg)
    dens = np.empty((lat_deg.size, 7))
    temp = np.empty((lat_deg.size, 2))
    for i in range(lat_deg.size):
        dens[i, :], temp[i, :] = venus_vts3.vts3(alt_km, 0, lat_deg[i], hour, f107a, f107, 48)

    var = 'lat_deg'
    atmo = xarray.Dataset(coords={var: lat_deg})
    atmo['CO2'] = (var, dens[:, 1])
    atmo['O'] = (var, dens[:, 2])
    atmo['CO'] = (var, dens[:, 3])
    atmo['He'] = (var, dens[:, 4])
    atmo['N'] = (var, dens[:, 5])
    atmo['N2'] = (var, dens[:, 6])

    atmo['Texo'] = (var, temp[:, 0])
    atmo['Tn'] = (var, temp[:, 1])

    atmo.attrs['hour'] = hour
    atmo.attrs['f107a'] = f107a
    atmo.attrs['f107'] = f107
    atmo.attrs['alt_km'] = alt_km

    return atmo
