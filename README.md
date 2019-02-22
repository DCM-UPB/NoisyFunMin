[![Build Status](https://travis-ci.com/DCM-UPB/NoisyFunMin.svg?branch=master)](https://travis-ci.com/DCM-UPB/NoisyFunMin)
[![Coverage Status](https://coveralls.io/repos/github/DCM-UPB/NoisyFunMin/badge.svg?branch=master)](https://coveralls.io/github/DCM-UPB/NoisyFunMin?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/9b819da118734bf9bdc80ba1c8cabf8b)](https://www.codacy.com/app/NNVMC/NoisyFunMin?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DCM-UPB/NoisyFunMin&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/dcm-upb/noisyfunmin/badge)](https://www.codefactor.io/repository/github/dcm-upb/noisyfunmin)

# NoisyFunMin

C++ Library for minimising noisy functions, such as integrals computed with the Monte Carlo method.

In `doc/` there is a user manual in pdf and a config for doxygen.

In `examples/` and `test/` there are examples and tests for the library.


Some subdirectories come with an own `README.md` file which provides further information.


# Supported Systems

Currently, we automatically test the library on Arch Linux (GCC 8) and MacOS (with clang as well as brewed GCC 8).
However, in principle any system with C++11 supporting compiler should work.


# Requirements

- CMake, to use our build process
- (optional) valgrind, to run `./run.sh` in `test/`
- (optional) pdflatex, to compile the tex file in `doc/`
- (optional) doxygen, to generate doxygen documentation in `doc/doxygen`


# Build the library

Copy the file `config_template.sh` to `config.sh`, edit it to your liking and then simply execute the command

   `./build.sh`

Note that we build out-of-tree, so the compiled library and executable files can be found in the directories under `./build/`.


# First steps

You may want to read `doc/user_manual.pdf` to get a quick overview of the libraries functionality. However, it is not guaranteed to be perfectly up-to-date and accurate. Therefore, the best way to get your own code started is by studying the examples in `examples/`. See `examples/README.md` for further guidance.
