#!/bin/sh
(
cd ../../build/examples
rm ./*.out
./ex3.exe
)
python plot.py
