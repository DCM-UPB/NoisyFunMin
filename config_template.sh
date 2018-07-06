#!/bin/bash

OS_NAME=$(uname)

# Library name
LIBNAME="nfm"

# C++ compiler
CC="g++"

# C++ flags (std=c++11 is necessary)
FLAGS="-std=c++11 -Wall -Werror"

# Optimization flags
OPTFLAGS="-O3 -flto"

# Debugging flags
DEBUGFLAGS="-g -O0"

# MCIntegrator++ Library
MCI_FOLDER="/...../MCIntegratorPlusPlus"
IMCI="-I${MCI_FOLDER}/src/"
LMCI="-L${MCI_FOLDER}"
LIBMCI="-lmci"

#FFNN Library (used in ex3)
FFNN_FOLDER="/...../FeedForwardNeuralNetwork"
IFFNN="-I${FFNN_FOLDER}/src/"
LFFNN="-L${FFNN_FOLDER}"
LIBNAMEFFNN="ffnn"
LIBFFNN="-lffnn"
