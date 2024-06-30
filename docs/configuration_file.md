# Configuration File

This configuration file is used to set various parameters for DC3DFWD. It must be
end with extension of `prm`. This configuration file is a simple text file where
each line represents a parameter for the program. The format of each line is as
follows:

```set <parameter_name> = <parameter_value>```

`set` is a keyword indicating that a parameter is being defined. `<parameter_name>`
is the name of the setting. It is case sensitive. `<parameter_value>` is the value
of the parameter. It can be a number, a string, or another type of value depending
on the specific setting. Blank lines and lines starting with a `#` are ignored.
Lines starting with a `#` are considered comments and can be used to add
explanations or notes to the file.

Below is a description of each parameter:

- `iprefix`: This sets the prefix for input files. The input files are the model
  files and the data file. Please refer to the [Model files](/docs/model_files.md)
  and [Data file](/docs/data_file.md) for details about the explanation of these files.

- `oprefix`: This sets the prefix for output files. The output files are the
  solution files and the response file. The solution files are `vtk` files
  containing the mesh, estimated error, and the solution. The response file is a
  text file containing the potential values at the observation points. The
  response file has extension of `rsp`. Please refer to the
  [Response file](/docs/response_file.md) for details about the response file.

- `fe_degree`: This sets the polynomial degree of finite elements used in the
  forward modeling. The default value is `1`.

- `solver`: This sets the type of preconditioner used to solve the linear
  system. The available preconditioners are `GMG` (Geometric Multigrid), `AMG`
  (Algebraic Multigrid). The default solver is `GMG`.

- `smoother_steps`: This sets the number of smoothing steps for the Geometric
  Multigrid solver. The default value is `1`.

- `max_iter`: This sets the maximum number of iterations for the iterative
  solver. The default value is `100`.

- `rtol`: This sets the relative tolerance for the iterative solver. The
  default value is `1E-8`.

- `n_adaptive_refinements`: This sets the number of adaptive mesh refinements to
  be performed in the forward modeling process. The default value is `0`.

- `refine_fraction`: This sets the fraction of cells to be refined in each
  adaptive mesh refinement. The default value is `0.1`, which means 10% of the
  cells will be refined in each adaptive mesh refinement step and will
  approximately double the number of cells.

- `coarsen_fraction`: This sets the fraction of cells to be coarsened in each
  adaptive mesh refinement. The default value is `0.0`.

- `n_global_refinements`: This sets the number of global mesh refinements to be
  performed prior to the adaptive mesh refinements. The default value is `0`.

Here is an example of a configuration file:

```text
set iprefix = halfspace
set oprefix = halfspace

set fe_degree = 1

set solver = GMG

set rtol = 1E-8
set max_iter = 100

set n_adaptive_refinements = 5
set n_global_refinements = 0
```
