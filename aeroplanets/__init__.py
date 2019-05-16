import logging
from pathlib import Path
import numpy as np
from typing import List, Tuple
import xarray
import xml.etree.ElementTree as etree
from datetime import datetime, timedelta

from . import venus  # noqa: F401


def filefind(path: Path) -> List[Path]:
    EXT = ('*.dat', '*.out')
    """find all data files in directory"""
    path = Path(path).expanduser()

    if path.is_file():
        return [path]

    flist = []
    for p in EXT:
        flist += list(path.glob(p))
# %% excise non-data files
    flist = [f for f in flist if f.name != 'emission_list.out']

    if not flist:
        raise FileNotFoundError(f'no files found in {path}')

    return flist


def loadaeroout(fn: Path, arr: xarray.Dataset) -> xarray.Dataset:
    # %% raw data
    data = np.loadtxt(fn, comments='#')
    if (data.ndim != 2 or np.allclose(data[:, 1:], 0.)):  # no data in file or all zero data
        logging.info(f'nothing in {fn}')
        return arr
# %% metadata
    lbl, cols = getmeta(fn, data.shape[1]-1)
# %% assemble
    """underlying coordinates are time, altitude, lat, lon"""
    lbl['filename'] = fn.name
    for i, c in enumerate(cols):
        if i == 0 and arr.alt_km.size == 0:
            arr['alt_km'] = data[:, 0]
        try:
            arr[c] = (('altkm',), data[:, i+1])
            # name=fn.name[:4]) # not stored in Dataset
            arr[c].attrs = lbl
        except ValueError:
            logging.error(f'{fn} shape inconsistent with other files, skipping.')

    return arr


def xmlparam(fn: Path) -> Tuple[datetime, dict]:
    """read data from simulation XML input file"""
    if not isinstance(fn, (str, Path)):
        raise ValueError('must specify simulation configuration .xml file')

    fn = Path(fn).expanduser()

    root = etree.ElementTree(file=fn).getroot()

    year = list(root.iterfind('planet/year'))[0].text
    md = list(root.iterfind('planet/day_of_year'))[0].text
    t = datetime.strptime(year+md, '%Y%j')  # year, month, day
    hod = float(list(root.iterfind('planet/UT'))[0].text)  # 0..24 UTC
    t += timedelta(hours=hod)

    glat = float(list(root.iterfind('planet/planet_lat'))[0].text)
    glon = float(list(root.iterfind('planet/planet_long'))[0].text)

    p = {'glat': glat, 'glon': glon}

    return t, p


def getmeta(fn: Path, ncol: int):
    fn = Path(fn).expanduser()
    """get metadata from raw file output"""
    lbl = {'filename': fn.name}
    cols = []
    with fn.open('r') as f:
        lbl['title'] = f.readline().strip()[2:]
        lbl['xlabel'] = f.readline().strip()[2:]
        lbl['ylabel'] = f.readline().strip()[2:]

        if fn.name.startswith('column'):  # TODO how to label wavelengths in these?
            for i in range(ncol):
                cols.append(fn.stem[6:]+'_'+str(i))
            return lbl, cols

        for i in range(ncol):
            line = f.readline().strip()
            if not line.startswith('#'):
                logging.warning(f'{fn} header may be read incorrectly')
                return lbl, cols
            cols.append(line[2:])

    return lbl, cols


def createds(xmlfn: Path) -> xarray.Dataset:
    """
    intializes the dataset, setting up the data coordinates
    """

    t, attrs = xmlparam(xmlfn)

    arr = xarray.Dataset(coords={'time': t,
                                 'alt_km': [],
                                 'lat': attrs['glat'],
                                 'lon': attrs['glon']})

    return arr


def convertdata(path: Path, xmlfn: Path, ofn: Path, types: tuple) -> xarray.Dataset:
    """raw text output to NetCDF4"""

    flist = filefind(path)

    data = createds(xmlfn)

    for f in flist:
        if not f.name[:4] in types:
            continue
        data = loadaeroout(f, data)

    print(len(data.variables), 'variables found in', flist[0].parent,
          'from', len(flist), 'files.')

    if ofn:
        ofn = Path(ofn).expanduser()
        if ofn.is_dir():
            raise OSError('must specify a NetCDF path/filename to write')
        print('writing', ofn)
        data.to_netcdf(ofn, mode='w')

    return data
