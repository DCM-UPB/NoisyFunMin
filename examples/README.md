# LEGEND OF THE EXAMPLES

Make sure the examples are compiled, by running `./build.sh` in the project root folder.
Execute an example by switching into one of the example folders and running `./run.sh`.
Note that the actual example executables reside inside the `build/examples/` folder under the project's root.

## Common

In the folder `common/` you will find a file named `ExampleFunctions.hpp` that contains some simple wave function
and Hamiltonian implementations that are shared between the examples.


## Example 1

`ex1/`: minimisation using ConjugateGradient of a quadratic function, without and with a noise.


## Example 2

`ex2/`: as `ex1`, but using the DynamicDescent method.


## Example 3

`ex3/`: as `ex1`, but using the Adam method.
