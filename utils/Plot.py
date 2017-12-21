#!/usr/bin/env python
from pathlib import Path
from matplotlib.pyplot import figure, draw,pause,show,close
import re
import numpy as np
from argparse import ArgumentParser


def filefind(path:Path):
    path = Path(path).expanduser()

    flist = []
    for p in ('*.dat','*.out'):
        flist += list(path.glob(p))

    return flist


def plotter(fn:Path, save:bool):
    """spawn new figure with data from file"""

    fig,ax,legende = _labeler(fn)
# %%
    data=np.loadtxt(fn, comments='#')
    if data.shape[1] < 2: # no data in file
        return

    ax.plot(data[:,1:],data[:,0])
#    for i in range(1, data.shape[1]):
#        ax.plot(data[:,i], data[:,0], label=legende[i-1])

    ax.set_xscale("log")
    ax.set_xlim(1e-2,None)
    ax.legend(legende[:-1],loc='best')
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

    flist = filefind(p.path)
    print('generating',len(flist),'plots.')
    for f in flist:
        try:
            plotter(f,p.save)
        except Exception as e:
            print(f,e)
        show()