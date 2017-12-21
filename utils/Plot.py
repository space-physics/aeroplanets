#!/usr/bin/env python
"""
Plot all output of aero1d program.
Python >= 3.5
"""

from pathlib import Path
import re
import numpy as np
from argparse import ArgumentParser

log = {'y':'Electron_precip.dat'}

def filefind(path:Path):
    path = Path(path).expanduser()

    if path.is_file():
        return [path]

    flist = []
    for p in ('*.dat','*.out'):
        flist += list(path.glob(p))

    flist = [f for f in flist if f.name != 'emission_list.out']

    if not flist:
        raise FileNotFoundError('no files found in {}'.format(path))

    return flist


def plotter(fn:Path, save:bool):
    """spawn new figure with data from file"""
    data=np.loadtxt(fn, comments='#')
    if data.ndim != 2: # no data in file
        return
    if np.allclose(data[:,1:], 0.): # all zero data
        return

    fig,ax,legende = _labeler(fn)
# %%

    ax.plot(data[:,1:],data[:,0])

    ax.set_xscale("log")
    if fn.name in log['y']:
        ax.set_yscale('log')

    ax.set_xlim(1e-2,None)
    ax.legend(legende[:data.shape[1]],loc='best')
    # %%
    if save:
        ofn = fn.with_suffix('.svg')
        fig.savefig(str(ofn))
        close(fig)


def _labeler(fn:Path):
    is_title=False
    is_x=False
    is_y=False
    legende=[]

    fig = figure()
    ax = fig.gca()

    with fn.open('r') as f:
        for lines in f:
            if(re.match("^#",lines)):
                if not is_title:
                    ax.set_title(lines[2:])
                    is_title=True
                else:
                    if not is_x:
                        ax.set_xlabel(lines[2:])
                        is_x=True
                    else:
                        if not is_y:
                            ax.set_ylabel(lines[2:])
                            is_y=True
                        else:
                            legende.append(lines[2:])

    return fig, ax, legende


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument("path", help="The aero1d output directory to plot from")
    p.add_argument("-save", help="write output image", action='store_true')
    p = p.parse_args()
# %% greatly speed up plotting when figures are just saved to disk
    if p.save:
        import matplotlib
        matplotlib.use('svg')

    from matplotlib.pyplot import figure,show,close
# %%
    flist = filefind(p.path)
    print('generating',len(flist),'plots.')
    for f in flist:
        try:
            plotter(f,p.save)
        except Exception as e:
            print(f,e)
        show()