import logging
from pathlib import Path
import numpy as np
import xarray as xr
import xml.etree.ElementTree as etree
from datetime import datetime,timedelta


def filefind(path:Path):
    EXT = ('*.dat','*.out')
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


def loadaeroout(fn:Path, arr:xr.Dataset):
# %% raw data
    data = np.loadtxt(fn, comments='#')
    if (data.ndim != 2 or np.allclose(data[:,1:], 0.)):  # no data in file or all zero data
        logging.info(f'nothing in {fn}')
        return arr
# %% metadata
    lbl,cols = getmeta(fn, data.shape[1]-1)
# %% assemble
    """underlying coordinates are time, altitude, lat, lon"""
    lbl['filename'] = fn.name
    for i,c in enumerate(cols):
        try:
            arr[c] = xr.DataArray(data[:,i+1],
                                  coords={'alt_km':data[:,0]},
                                  dims=['alt_km'],)
                                  #name=fn.name[:4]) # not stored in Dataset
            arr[c].attrs=lbl
        except ValueError:
            logging.error(f'{fn} shape inconsistent with other files, skipping.')

    return arr


def xmlparam(fn:Path):
    """read data from simulation XML input file"""
    fn = Path(fn).expanduser()

    root = etree.ElementTree(file=fn).getroot()

    year = list(root.iterfind('planet/year'))[0].text
    md = list(root.iterfind('planet/day_of_year'))[0].text
    t = datetime.strptime(year+md,'%Y%j') # year, month, day
    hod = float(list(root.iterfind('planet/UT'))[0].text) # 0..24 UTC
    t += timedelta(hours=hod)

    glat = float(list(root.iterfind('planet/planet_lat'))[0].text)
    glon = float(list(root.iterfind('planet/planet_long'))[0].text)

    p = {'glat':glat, 'glon':glon}

    return t,p


def getmeta(fn:Path, ncol:int):
    fn = Path(fn).expanduser()
    """get metadata from raw file output"""
    lbl = {'filename':fn.name}
    cols = []
    with fn.open('r') as f:
        lbl['title']  = f.readline().strip()[2:]
        lbl['xlabel'] = f.readline().strip()[2:]
        lbl['ylabel'] = f.readline().strip()[2:]

        if fn.name.startswith('column'): # TODO how to label wavelengths in these?
            for i in range(ncol):
                cols.append(fn.stem[6:]+'_'+str(i))
            return lbl,cols

        for i in range(ncol):
            line = f.readline().strip()
            if not line.startswith('#'):
                logging.warning(f'{fn} header may be read incorrectly')
                return lbl, cols
            cols.append(line[2:])

    return lbl,cols


def createds(xmlfn:Path, flist:list):

    data = np.loadtxt(flist[0], comments='#')

    t, attrs = xmlparam(xmlfn)

    arr = xr.Dataset(coords = {'time': t,
                                'alt_km': data[:,0],
                                'lat':attrs['glat'],
                                'lon':attrs['glon']})

    arr['filelist'] = [f.name for f in flist]

    return arr


def convertdata(fn:Path, xmlfn:Path, ofn:Path, types:tuple):
    """raw text output to NetCDF4"""

    flist = filefind(fn)

    data = createds(xmlfn, flist)

    for f in flist:
        if not f.name[:4] in types:
            continue
        data = loadaeroout(f, data)

    print(len(data.variables), 'variables found in', flist[0].parent,
          'from', len(flist), 'files.')

    if ofn:
        ofn = Path(ofn).expanduser()
        print('writing',ofn)
        data.to_netcdf(str(ofn), mode='w')


    return data