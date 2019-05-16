from pathlib import Path
from matplotlib.pyplot import figure, close
#
log = {'y': ('Electron_precip.dat',),
       'x': ()}


def plotalt(ax, data, d):
    ax.plot(data[d], data.alt_km, label=d)
    ax.set_title(data[d].attrs['title'])
    ax.set_xlabel(data[d].attrs['xlabel'])
    ax.set_ylabel(data[d].attrs['ylabel'])
    ax.grid(True)
    ax.legend(loc='best')


def plotter(data, typ: str, save: Path):
    """spawn new figure with grouped data from file"""

    fig = figure(figsize=(18, 12))
    ax = fig.gca()
# %% plots
    attrs = None
    for d in data.data_vars:
        if data[d].attrs['filename'][:4] != typ:
            continue

        if d.startswith('Electron (TEC'):
            plotalt(figure(figsize=(18, 12)).gca(), data, d)
            continue

        ax.plot(data[d], data.alt_km, label=d)
        attrs = data[d].attrs  # use last match

    if attrs is None:
        print('no', typ, 'found')
        close(fig)
        return

    if not typ == 'iono':
        ax.set_xscale("log")
        ax.set_xlim(1e-2, None)
    ax.grid(True)

    ax.set_title(attrs['title'])
    ax.set_xlabel(attrs['xlabel'])
    ax.set_ylabel(attrs['ylabel'])

    ax.legend(loc='best')
    # %%
    if save:
        save = Path(save).expanduser()
        ofn = Path(save / (typ + '.png'))
        print('writing', ofn)
        fig.savefig(ofn)
        close(fig)
