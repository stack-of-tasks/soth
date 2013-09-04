soth
====

[![Build Status](https://travis-ci.org/stack-of-tasks/soth.png)](https://travis-ci.org/stack-of-tasks/soth)
[![Coverage Status](https://coveralls.io/repos/stack-of-tasks/soth/badge.png?branch=master)](https://coveralls.io/r/stack-of-tasks/soth?branch=master)


Setup
-----

To compile this package, it is recommended to create a separate build
directory:

    mkdir _build
    cd _build
    cmake [OPTIONS] ..
    make install

Please note that CMake produces a `CMakeCache.txt` file which should
be deleted to reconfigure a package from scratch.


### Dependencies

The matrix abstract layer depends on several packages which
have to be available on your machine.

 - Libraries:
   - eigen3
 - System tools:
   - CMake (>=2.6)
   - pkg-config
   - usual compilation tools (GCC/G++, make, etc.)
