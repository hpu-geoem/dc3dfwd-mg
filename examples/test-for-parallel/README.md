# Steps to run the example

This example demonstrates how to run the forward modeling process the half-space
model with different number of processes.

Modify the multigrid method and stretch factor in halfspace.prm and mkmdl.py, respectively.

Step into each subdirectory and run the following command to generate the model
and data files, and run the forward modeling process:

```shell
python mkmdl.py
mpirun -np <number of processors> ../../dc3dfwd halfsapce.prm
```

Please replace `<number of processors>` with the number of processors you want to use.
