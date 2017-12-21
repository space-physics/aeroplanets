# Trans++: C++ implementation of the Trans* programs.


## Prereqs

* Mac: `brew install make gcc boost doxygen`
* Linux: `apt install make g++ libboost-filesystem-dev libboost-regex-dev libboost-thread-dev libboost-program-options-dev doxygen graphviz`

Plots use Python and Matplotlib as is customary.

## Build
```sh
./prepare_conf.sh

./configure --enable-maintainer-mode --enable-debug

make

make check

make install
``` 

may give some errors on the newest version of boost.


## Usage

Aurora example

1. run sim
   ```sh
   cd data/Earth 

   ./aero1d AuroraEarthFairbanks.xml
   ```
2. Plot
   ```sh
   ./SortieAurora/bin/plotte`
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







