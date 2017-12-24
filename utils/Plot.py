#!/usr/bin/env python
"""
Convert and/or Plot all output of aero1d program.

./utils/Plot.py data/Earth/SortieAurora/ data/Earth/AuroraEarthFairbanks.xml test.nc

./utils/Plot.py test.nc
"""
from pathlib import Path
import xarray as xr
from matplotlib.pyplot import show
import seaborn as sns
sns.set_context('talk',font_scale=1.2)
#
from aeroplanets import convertdata
from aeroplanets.plots import plotter
#
ALONGB = ['chem','colu','extr','iono','neut','prod'] # MUST be list, not tuple

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("path", help="The aero1d output directory to plot from (or NetCDF4 file to load)")
    p.add_argument('xmlfn', help='path to simulation XML config file',nargs='?')
    p.add_argument('ofn', help='NetCDF4 file to write',nargs='?')
    p.add_argument("-saveplots", help="write figures", action='store_true')
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