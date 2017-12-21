#!/usr/bin/env python
from pathlib import Path
from matplotlib.pyplot import figure, draw,pause,show
import re
import numpy as np
from argparse import ArgumentParser


def filefind(path):
    path = Path(path).expanduser()

    flist = []
    for p in ('*.dat','*.out'):
        flist += list(path.glob(p))

    return flist


def plotter(fn):
    """spawn new figure with data from file"""


    fig,ax,legende = _labeler(fn)
# %%
    data=np.loadtxt(fn)

    for i in range(1, len(data[0])):
        ax.plot(data[:,i], data[:,0], label=legende[i-1])

    ax.set_xscale("log")
    ax.set_xlim(1e-2,None)
    ax.legend(loc='best')
    # %%
    if not p.outfn:
        imagename = p.path+".png"

    fig.savefig(imagename)


def _labeler(fn):
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
    p.add_argument("-o","--outfn", help="The output image")
    p = p.parse_args()

    flist = filefind(p.path)
    print('generating',len(flist),'plots.')
    for f in flist:

        try:
            plotter(f)
        except Exception as e:
            print(f,e)
        draw()
        pause(0.01)
    show()