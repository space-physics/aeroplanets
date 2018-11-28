# Trans++: C++ implementation of the Trans* programs.

is from [Guillaume Gronoff](https://scholar.google.com/citations?user=e2RfvmYAAAAJ) of NASA, a former student of P. Blelly, Transcar author.
Michael Hirsch changed the broken autoconf build system to Cmake and redid the plotting utilities.



It builds easily with Cmake. 
It can be compared with other models like 
[Glow](https://www.github.com/scivision/glowaurora) or 
[Transcar](https://www.github.com/scivision/transcar).

The simulation is configured completely with XML.
The runtime is similar to Transcar (3-5 minutes on a laptop). 
The output is text files, for which I have fast Python parsers that could dump to HDF5 etc. as desired. 
Over 70 output files are created per simulation, covering production rates for numerous neutral and ion species, column emissions from UV to IR, energy grid from 1 eV ~ 30 keV, effects of secondary electroncs are considered.  

Inputs include:
* photon flux (photoionization)
* proton flux (proton impact ionization) 
* electron flux (electron impact ioniziation) -- parameterization or read your e- flux file.



## Prereqs

* Mac: `brew install make cmake gcc boost doxygen openblas`
* Linux: `apt install make cmake g++ gfortran libopenblas-dev libboost-filesystem-dev libboost-regex-dev doxygen graphviz`
* BSD: `pkg install make cmake gcc boost-libs doxygen openblas`
* Windows: suggest [Windows Subsystem for Linux](https://www.scivision.co/install-windows-subsystem-for-linux/)

Plotting uses Python &ge; 3.5 and Matplotlib.

## Build
```sh
cd bin
cmake ../src
make -j -l2  # don't omit the -l2 or the computer will crash on compile due to excess resource use
```

## Usage
Currently there is a bug where if the output directory doesn't exist, the problem just keeps computing...but doesn't write any output at all.

Aurora example

1. run sim
   ```sh
   cd ../data/Earth/
   ../../bin/aero1d AuroraEarthFairbanks.xml
   ```
   the output appears in `data/Earth/SortieAurora`
2. Plot--automatically iterates over all files in output directory.
   ```sh
   ./utils/Plot.py data/Earth/SortieAurora/ -save
   ```


`SortieAuroraCompar` is reference data to check if the compilation/computation was successful.

### Precipitation
Alfvenic aurora and other structured aurora users may be interested in configuring electron precipitation flux characteristics.
This is configured in the input XML file as follows:

1. `<use_precipitation />` must exist to enable electron precipitation, and the parameters are enclosed by `<precipitation>` `</precipitation>` XML tags.
2. Electron precipitation must be read from a file or computed by a simple parameterization.

`use_model`  (created in `math/mathflux.cpp`)
* 0: Null (none)
* 1: Maxwellian
* 2: Gaussian
* 3: Dirac (pseudo-monoenergetic)


* `E0`     Characteristic energy [eV]
* `entot`  Total flux [ergs]
* `isotro` 0: non-isotropic.  1: isotropic.
* `powlaw` power-law shaping 0 or 1


### Plotting
Several programs related to concurrent running of large numbers of simulations are in `utils/`, for advanced users only.
Most users will simply use `Plot.py` to make line plots of text output files, using headers.

## Documentation

1. compile docs
   ```sh
   doxygen Doxyfile
   ```
2. web browser: `doc/html/index.html`
3. To create `/doc/latex/refman.pdf`:
   ```sh
   cd doc/latex
   make
   ```

