# Steps to run the example

This example demonstrates how to run the forward modeling process the half-space
model with different orders of elements and different refinement method using GMG.

Step into each subdirectory and run the following command to generate the model
and data files, and run the forward modeling process:

```shell
python mkmdl.py
mpirun -np <number of processors> ../../dc3dfwd GMG-adaptive-refinement-1.prm
mpirun -np <number of processors> ../../dc3dfwd GMG-adaptive-refinement-3.prm
mpirun -np <number of processors> ../../dc3dfwd GMG-global-refinement-1.prm
mpirun -np <number of processors> ../../dc3dfwd GMG-global-refinement-3.prm
```

Please replace `<number of processors>` with the number of processors you want to use.