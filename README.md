# Trans++: C++ implementation of the Trans* programs.


## Prereqs

* Mac: `brew install make cmake gcc boost doxygen openblas`
* Linux: `apt install make cmake g++ libopenblas-dev libboost-filesystem-dev libboost-regex-dev libboost-thread-dev libboost-program-options-dev libboost-test-dev doxygen graphviz`

Plots use Python and Matplotlib as is customary.

## Build
```sh
cd bin
cmake ../src
make -j4
```

### Build (old way, unstable)
I don't use this anymore as it was unreliable and not forward compatible.

```sh
./prepare_conf.sh

./configure --enable-maintainer-mode --enable-debug

make

make check

make install
```

may give some errors on the newest version of boost.


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
3. `SortieAuroraCompar` is reference data to check if the compilation/computation was successful.

### Precipitation
Alfvenic aurora and other structured aurora users may be interested in configurings electron precipitation flux characteristics.
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
Several plotting programs exist in `utils`.

* Plot.py: line plots of text output files, using headers.

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







