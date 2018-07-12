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