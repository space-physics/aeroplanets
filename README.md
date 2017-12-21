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
```
Some usual configurations are listed in  ./DifferentConfigurations

#>make check 
may give some errors on the newest version of boost.


When the program is compiled and installed. You can go to data/Earth 
and run 
#>aero1d AuroraEarthFairbanks.xml
Then, you can go to SortieAurora and use #> plotte (if copied in the bin directory) to plot the different outputs.
SortieAuroraCompar is pre-computed and allows to check if the compilation/computation was successful.


## Documentation
```sh
doxygen Doxyfile
```

`doc/` must have `html/` and `latex/` directory.

To browse the configuration with an internet browser, open
`doc/html/index.html`

To have a pdf document, go to the latex directory and type make.
If you have the good latex libraries, it will create a refman manual.








