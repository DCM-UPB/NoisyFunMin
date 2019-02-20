[![Build Status](https://travis-ci.com/DCM-UPB/NoisyFunMin.svg?branch=master)](https://travis-ci.com/DCM-UPB/NoisyFunMin)
[![Coverage Status](https://coveralls.io/repos/github/DCM-UPB/NoisyFunMin/badge.svg?branch=master)](https://coveralls.io/github/DCM-UPB/NoisyFunMin?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/9b819da118734bf9bdc80ba1c8cabf8b)](https://www.codacy.com/app/NNVMC/NoisyFunMin?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DCM-UPB/NoisyFunMin&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/dcm-upb/noisyfunmin/badge)](https://www.codefactor.io/repository/github/dcm-upb/noisyfunmin)

# NoisyFunMin

C++ Library for minimising noisy functions, such as integrals computed with the Monte Carlo method.

In `doc/` there is a user manual in pdf.

In `examples/` there are several examples.



# Build the library

We use the CMake build system, so you need to have it on your system to build the library out of the box.
Then copy the file `config_template.sh` to `config.sh`, edit it to your liking and then simply execute the command

   `./build.sh`

Note that we build out-of-tree, so the compiled library and executable files can be found in the directories under `./build/`.
