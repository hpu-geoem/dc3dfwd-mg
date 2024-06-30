# DC3DFWD

## About

DC3DFWD is a parallel C++ program designed for the forward modeling of the 3D
direct current (DC) resistivity method using the adaptive finite element method
and the geometric multigrid method. It is designed to run on high
performance computing (HPC) systems. The main features of DC3DFWD including:

- Implementation of high-order adaptive finite elements to adaptively refine the
  mesh according to the go-oriented error estimator
- Efficient geometric multigrid solver on locally refined meshes
- Parallel computation using MPI
- Discretization of complex structures using octree meshes

## Installing Prerequisites and Building

DC3DFWD is developed based on the open-source library
[deal.II](https://www.dealii.org/), which is a C++ program library that
provides building blocks for the finite element method.
We provide two ways to install deal.II and its dependencies: using `Docker` and
using `spack`. Among the two approaches, we recommend using `Docker`, since it
is the most robust and portable way to build DC3DFWD.

Note that DC3DFWD has only been tested on Linux and macOS systems. For
Windows users, we recommend using Docker, WSL1/WSL2 or a virtual machine.

Please extract the source code of DC3DFWD into a directory
`/path/to/dc3dfwd-sources` before installing the prerequisites and building.
In the following sections, we will use `/path/to/dc3dfwd-sources` to refer to
the source directory of DC3DFWD.

### Using Docker

If the user is familiar with `Docker`, we provide a Docker image that contains
all the dependencies of DC3DFWD. For macOS and Windows users, please install
[Docker Desktop](https://www.docker.com/products/docker-desktop) first. For
Linux users, please install `Docker Engine` following the
[installation instructions](https://docs.docker.com/engine/install).

There is an alternative Docker client called `OrbStack` for macOS users.
It is claimed to be faster than `Docker Desktop` and has a better user
interface. Please refer to [OrbStack](https://orbstack.dev) for more details.

After installing `Docker`, the user can pull the Docker image from
[Docker Hub](https://hub.docker.com/r/adamqc/dc3dfwd-deps) using the following
command:

```shell
docker pull adamqc/dc3dfwd-deps:latest
```

Then run the Docker container and mount the source directory of DC3DFWD to
`/mnt` using the following command:

```shell
docker run -it --rm -v /path/to/dc3dfwd-sources:/mnt adamqc/dc3dfwd-deps:latest
```

The above command will start a shell inside the Docker container. Now the user
can build DC3DFWD:

```shell
cd /mnt && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

### Using spack

To build deal.II and its dependencies from source, we recommend using
[spack](https://spack.io), which is a package manager for supercomputers,
Linux, and macOS. Please refer to the [spack documentation](https://spack.readthedocs.io/en/latest/)
for more details.

Before we install `spack`, we need to install a C++ compiler and other
dependencies. For debian based Linux distributions, e.g., Ubuntu, Mint, etc.,
please use the following command to install GCC and other necessary packages:

```shell
sudo apt install -y build-essential gfortran git python3
```

For RedHat based Linux distributions, e.g., RHEL, Centos, etc., please use the
following command:

```shell
sudo yum install -y gcc gcc-c++ gcc-gfortran git python3
```

On macOS, the user can install `Xcode` from the App Store, or install
`Command Line Tools for Xcode` manually.

Note that `deal.II` requires a C++ compiler that supports C++17 standard, e.g.,
GCC >= 9.0, Clang >= 10.0, and AppleClang >= 12.0. The user can check the
version of the compiler using `g++ --version` or `clang++ --version`.
Also note that AppleClang >= 15.0 has a problem compiling `deal.II`, so we
recommend using `Xcode` or `Command Line Tools for Xcode` prior to version 14.3.

The next step is to install `spack`. We need the developing version of spack
since it contains the latest version of these packages. The following commands
can be used to install `spack` and `deal.II`:

```shell
# Clone the spack repository
git clone https://github.com/spack/spack.git
# Initialize spack environment
source ./spack/share/spack/setup-env.sh
# Install deal.II and its dependencies
spack install dealii@master+mpi+petsc+python~p4est~arborx~arpack~slepc~gsl~adol-c~hdf5~gmsh~sundials~oce~cgal~assimp~symengine~examples~ginkgo~threads~muparser~vtk build_type=Release ^python@3.10
```

Note that the `+python` and `+petsc` options are required since DC3DFWD uses
Python to generate the mesh and PETSc to solve the linear system. This command
may take a while to finish since it will build all the packages from source,
please be patient.

Once all the dependencies are installed, the following commands can be used to
build DC3DFWD:

```shell
# Go to the source directory of DC3DFWD
cd /path/to/dc3dfwd-sources
# Load environment variables of deal.II
spack load dealii
# Create a build directory
mkdir build
# Go to the build directory
cd build
# Configure
cmake -DCMAKE_BUILD_TYPE=Release ..
# Build
make
```

The above commands will generate an executable file `dc3dfwd` in the source
directory, which is the main program of DC3DFWD.

## Usage

To run the forward modeling process, the user needs to provide a configuration file
to specify the parameters. Then use the flowing command to run DC3DFWD:

```shell
mpirun -np <number of processors> /path/to/dc3dfwd path/to/options.prm
```

For details about the options file, please refer to the [Configuration file](/docs/configuration_file.md).

## License

DC3DFWD is distributed under the MIT License.

## Contributing

Users are encouraged to open an issue for questions or bugs. Pull requests for
any enhancements are also welcome.
