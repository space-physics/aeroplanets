#!/usr/bin/env python
from aeroplanets.plots import plotter
from aeroplanets import convertdata
"""
Convert and/or Plot all output of aero1d program.

python utils/Plot.py build/SortieAurora/ data/Earth/AuroraEarthFairbanks.xml test.nc

python utils/Plot.py test.nc
"""
from pathlib import Path
import xarray as xr
from matplotlib.pyplot import show
try:
    import seaborn as sns
    sns.set_context('talk', font_scale=1.2)
except ImportError:
    print('recommend "conda install seaborn" for better looking plots')
#
#
ALONGB = ['chem', 'colu', 'extr', 'iono', 'neut', 'prod']  # MUST be list, not tuple

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("path", help="The aero1d output directory to plot from (or NetCDF4 file to load)")
    p.add_argument('xmlfn', help='path to simulation XML config file', nargs='?')
    p.add_argument('ofn', help='NetCDF4 file to write', nargs='?')
    p.add_argument("-p","--saveplots", help="path to save plots to")
    p = p.parse_args()
# %% along-B data
    path = Path(p.path).expanduser()
    if path.suffix == '.nc':
        data = xr.open_dataset(path)
    else:
        data = convertdata(path, p.xmlfn, p.ofn, ALONGB)
# %% plot
    for q in ALONGB:
        plotter(data, q, p.saveplots)

    show()
