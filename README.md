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
   ./aero1d AuroraEarthFairbanks.xml 
   ```
   the output appears in `data/Earth/SortieAurora`
2. Plot
   ```sh
   ./utils/Plot.py data/Earth/SortieAurora
   ```
3. `SortieAuroraCompar` is pre-computed and allows to check if the compilation/computation was successful.


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







