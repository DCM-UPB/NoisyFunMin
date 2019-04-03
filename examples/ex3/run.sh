#!/bin/sh
(
cd ../../build/examples
rm *.out
ln -sf ../../examples/ex3/plot.py
./ex3.exe
)
python plot.py
