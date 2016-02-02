Copyright 2012-2016 Mario Schröck, Hannes Vogt (culgt.com)

cuLGT comes with stand-alone lattice gaugefixing applications in lib/gaugefixing and
an (very experimental!) interface to the MILC framework in lib/milcinterface.


# SUPPORTED GAUGEFIXING #
This version (cuLGT 2.0) supports gaugefixing on a single GPU for

* **Landau gauge**, **Coulomb gauge** 
* single (SP), double (DP) and mixed precision (MP)
* SU(2) and SU(3)

If you are interested in applications for **maximally Abelian gauge** or **direct maximal center gauge**,
please contact one of the authors.

For multi-GPU support (Landau gauge) please use the culgt1 branch in the github repository.

# USING THE STANDALONE APPLICATIONS #
The standalone applications are located in lib/gaugefixing:

### DEPENDENCIES ###
* boost::program_options
* Random123 (is included in cuLGT, see license in include/external/Random123): https://www.deshawresearch.com/resources_random123.html
* google test and google mock [only for unit tests]
* c-lime and tinyxml [for ILDG file support (see below)]

### COMPILER/LIBRARY SUPPORT ###
Tested with

* NVCC version 6.5, 7.0, 7.5
* host compiler g++ version 4.8.3, 4.9.0
* BOOST version 1.50, 1.53, 1.58, 1.59, 1.60 
* -std=c++11 (and without)

Not all combinations work:

* BOOST 1.60 has a problem with NVCC if std != c++11 (http://stackoverflow.com/q/34959032/5085250)
* NVCC 6.5 does not allow g++-4.9 as host compiler

Some combinations will emit more warnings than others.

### BUILD ###

cuLGT comes with CMake build files. All relevant parameters can be set using CMake variables.

Mandatory parameters are (they are initialized to a default value, but most probably they won't fit):

* `CUDA_ARCH`: set to sm_XX according to your hardware
* `SUN`: set to 2 for SU(2) or 3 for SU(3)
* `CULGT_HOME`: path to the cuLGT root directory (alternative: set CULGT_HOME as environmental variable)

Optional parameters:

* `DISABLE_AUTOTUNE` (on/OFF): the optimal kernel configuration is tuned at runtime; the correspondant
   kernels are generated at compile-time which slows down compilation; for debugging purposes one might
   limit the compilation to one kernel with this option.
* `CULGT_ALLOW_C++11` (ON/off): possibility to disable c++11 compilation (automatically disabled for older compilers)
* `USE_TIMESLICEPATTERN` (on/OFF): run Landau gauge in timeslices

If you are not familiar with cmake, just run the following commands from `lib/gaugefixing`:


```
mkdir -p build
cd build
cmake .. -DCUDA_ARCHITECTURE=<your arch> -DSUN=<2 or 3> -DCULGT_HOME=<path to cuLGT root>
make
```

 
### RUN ###

To test the application on a random lattice run

```
./LandauGaugeSP --sethot 1
```

The application uses boost::program_options to manage command line options (or configuration file options).
Type ./LandauGaugeSP --help for a list of all options. Here we list the most important options:

* `--nx --nt (--ny --nz)`: lattice dimensions; if ny, nz is ommitted they default to nx)
* `--filetype`: MDP, NERSC, ILDG file types are supported
* `--fbasename`: path to a common name of all configurations, i.e. the part before numbering starts
* `--fextension`: file extension (if no extension is given it will default to a commonly used extension for the filetype)
* `--fnumberformat`: # of (min.) digits to use in the numbering
* `--reinterpret` FLOAT/DOUBLE: if you want to use the DP version of the program for a SP configuration you need to give FLOAT,
   i.e. the precision of the configuration, and vice versa.
* `--seed` for RNG (used in simulated annealing)

Example:
The configurations are called `/myconfigs/config_n16t16_000.dat`, filetype is MDP in DP (fermiqcd.net)

```
./LandauGaugeSP --nx 16 --nt 16 --fbasename /myconfigs/config_n16t16_ --fnumberformat 3 --filetype MDP --fextension .dat --reinterpret DOUBLE
```

or

```
./LandauGaugeDP --nx 16 --nt 16 --fbasename /myconfigs/config_n16t16_ --fnumberformat 3 --filetype MDP --fextension .dat
```

# TESTS #

The CMake interface for the unit tests needs some cleanup...
For the tests you need to install google test and google mock (https://github.com/google/googletest).
Then execute (in `test/build/`)
`cmake .. -DGTEST_HOME=<path to gtest home> -DGMOCK_HOME=<path to gmock home>`

The CMake script is looking for the static libraries `libgmock.a` in `${GMOCK_HOME}/build/` and `libgtest.a` in `${GTEST_HOME}/build/libgtest.a.`
If you use other locations specify `GTEST_LIBDIR` and `GMOCK_LIBDIR`.

After compilation you will find 3 executables: `test_all` (all tests), `test_cuda` (all cuda tests) and `test_host` (all host tests)

Additionally, you can compile the unit tests in each subdirectory separately (to save compilation time).

If you encounter any problems, please contact one of the authors.


# FILETYPES #

Some filetypes are not fully supported.
Fully supported:

* MDP
* HIREP

Partly supported:

* NERSC: If a NERSC file is loaded, the header information is used to write the file
          If a NERSC file is written (from a different file), only parts of the header information is written: DATATYPE, DIMENSION_X, FLOATING_POINT; other
          header information is set to 0 or n/a.
* ILDG:  Needs c-lime (https://github.com/usqcd-software/c-lime) and tinyxml (should be provided by repositories of your OS) otherwise the support is disabled
          at compile-time.
          If a file is written without prior load, the file will be damaged because the xml information is not written.

The limitations are only important if you choose to have a different output filetype compared to your input configuration.

For other (commonly used) file types we recommend the QCDUTILS by Massimo Di Pierro, see http://arxiv.org/abs/1202.4813 and https://code.google.com/archive/p/qcdutils/


# MILC INTERFACE #

If you want to use the MILC interface, please contact one of the authors.


# Known issues #

* (Very) large lattices are not supported at this stage for compute capability < 3.5. As a workaround you may disable the use of textures. Feel free to contact one of the authors.
* SP and MP use the same tune cache. If you are switching between SP and MP consider using separate directories for the application.