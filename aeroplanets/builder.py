from pathlib import Path
import os
import numpy.f2py
import importlib


def build(srcfn: Path, modname: str):
    """
    build Fortran module
    """
    srcfn = Path(srcfn).expanduser()
    if not srcfn.is_file():  # necessary check on Windows
        raise FileNotFoundError(f'{srcfn} is not a file')

    src = srcfn.read_text()
    opts = ['--quiet', '--f77flags="-w"']
    if os.name == 'nt':
        opts += ['--compiler=mingw32']

    numpy.f2py.compile(src, modname, opts)

    # verify that module was built and is importable
    importlib.import_module(modname)
