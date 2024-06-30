# Steps to run the example

This example demonstrates how to run the forward modeling process the 3-layers 
model using AMG and GMG.

Step into each subdirectory and run the following command to generate the model
and data files, and run the forward modeling process:

```shell
python mkmdl.py
mpirun -np <number of processors> ../../dc3dfwd AMG-global-refinement.prm
mpirun -np <number of processors> ../../dc3dfwd GMG-global-refinement.prm
```

Please replace `<number of processors>` with the number of processors you want to use.