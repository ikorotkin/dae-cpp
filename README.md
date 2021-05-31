# dae-cpp

[![Build Status](https://travis-ci.com/ikorotkin/dae-cpp.svg?branch=master)](https://travis-ci.com/ikorotkin/dae-cpp)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4aa33eb3a2834808a6cd1b81e0d8cc23)](https://www.codacy.com/app/ikorotkin/dae-cpp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ikorotkin/dae-cpp&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3256507.svg)](https://doi.org/10.5281/zenodo.3256507)

A simple but powerful C++ solver for Differential Algebraic Equation (DAE) systems.

## What is dae-cpp

A cross-platform, parallel C++ library for solving user-defined, stiff systems of DAEs (an initial value problem). The system may contain both differential and algebraic equations and can be written in the following matrix-vector form:

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathbf%7BM%7D%5Cfrac%7Bd%5Cmathbf%7Bx%7D%7D%7Bdt%7D%3D%5Cmathbf%7Bf%7D%5Cleft%20%28%20%5Cmathbf%7Bx%7D%20%5Cright%20%29">
</p>

where Mass matrix **M** can be singular, and the RHS **f**(**x**) is a nonlinear function of a real vector **x** and time *t*.

For the numerical integration the solver uses implicit [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) (Backward Differentiation Formula) method of orders 1-6 (can be defined by a user) with adaptive time stepping.

### How does it work

BDF time stepper reduces the original DAE system to a system of nonlinear equations that is solved using iterative [Newton root-finding algorithm](https://en.wikipedia.org/wiki/Newton%27s_method). Each Newton iteration a system of linear algebraic equations is solved using Parallel Direct Sparse Solver ([Intel MKL PARDISO](https://software.intel.com/en-us/mkl-developer-reference-c-intel-mkl-pardiso-parallel-direct-sparse-solver-interface)). The sparse solver performs 3 steps: reordering and symbolic factorization of Jacobian matrix, then numerical factorization, and then back substitution + iterative refinement. Finally, depending on the convergence rate of the Newton method, variability of the solution and user-defined accuracy, the DAE solver may adjust the time step and initiate a new iteration in time.

### The main features of the solver

-   Can resolve DAE systems of 10<sup>8</sup> equations and even more (depending on the Jacobian matrix sparsity and machine's RAM).
-   A user can provide analytical Jacobian matrix for better performance or use built-in parallel function provided by the solver to estimate numerical Jacobian.
-   Utilises all available cores on the machine for better performance (this can be overridden by a user).
-   Allows a user to adjust most of the parameters related to the solution process in order to achieve better accuracy and performance.
-   A user can get access to the solution at each time step by overriding Observer function (this is optional).
-   The library provides a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module for plotting.
-   The user-defined RHS, Mass matrix and Jacobian can be saved to a file for debugging or visualisation if needed.
-   Easy-to-follow examples (see, for example, [simple_dae.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/simple_dae/simple_dae.cpp), [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp) or [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp)) to kick-start the user's project.

## Installation

This is a cross-platform software that works on Linux (e.g. Ubuntu), Windows and macOS. The main library (the DAE solver itself) and all examples have only one external dependency: [Intel Math Kernel Library](https://software.intel.com/en-us/mkl), a fast and very well optimised math library. So the first step in the installation process is to download and install Intel MKL: [Linux](https://software.intel.com/en-us/mkl/choose-download/linux), [Windows](https://software.intel.com/en-us/mkl/choose-download/windows), [macOS](https://software.intel.com/en-us/mkl/choose-download/macos). Note that you may need to register in order to download the library. When asked, choose Intel Math Kernel Library for your OS and the Full Package.

An alternative and probably the most convenient way to download and install Intel MKL on Ubuntu (using APT Repository) is the following.

Install the GPG key for the repository:

```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```

Add the APT Repository:

```bash
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
```

Update the list of packages and install the library:

```bash
sudo apt-get update
sudo apt-get install intel-mkl-2019.5-075
```

This will install Intel MKL 2019.5. The list of all available versions and more information can be found [here](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo).

Note the latest versions of Intel MKL (2020) may produce a lot of run-time warnings. This is a known issue, the only workaround is to suppress them exporting the following variable:

```bash
# Suppress MKL run-time warnings (to fix a known issue of MKL 2020)
export KMP_WARNINGS=0
```

### Linux

On Linux make sure you have `git`, `cmake` and `g++` installed:

```bash
sudo apt-get install g++ cmake cmake-curses-gui git
```

In order to enable plotting (optional), `python3`, `matplotlib` and `numpy` should be installed:

```bash
sudo apt-get install python3 python3-dev python3-numpy python3-matplotlib
```

Then download dae-cpp library:

```bash
git clone https://github.com/ikorotkin/dae-cpp.git
```

The easiest way to install the library and compile all examples is just to create the build directory, then execute `cmake` (providing installation path) and `make`:

```bash
cd dae-cpp
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/install_path
make -j4
make install
```

where `/install_path` is the user-defined path where the package should be installed.

Note that `cmake` will try to find Intel MKL at its default location `/opt/intel/mkl` or according to `MKLROOT` environment variable. If the installation path is different, please provide MKL root path with the following `cmake` option: `-DDAE_MKL_DIR=/path_to_intel_mkl`.

Instead of `cmake -DCMAKE_INSTALL_PREFIX=/install_path ..` you might consider using `ccmake ..`, a GUI for `cmake` that will allow you to see all the options available before building the solver.

#### Test the solver

The DAE solver can perform a quick self test. To build the test, dae-cpp should be installed with `DAE_TEST=ON` option (it is ON by default). To start the test, from the build directory execute `ctest`:

```bash
ctest
```

During this test the solver will solve DAE systems from [examples](https://github.com/ikorotkin/dae-cpp/tree/master/examples) directory using analytical (if available) and numerical Jacobians, and then compare the results with the reference solutions.

#### More building options

-   `DAE_LONG_INT` - Use long integer representation for huge systems (more than ~10<sup>7</sup> equations). This option is OFF by default. For relatively small systems it is recommended to leave it OFF.
-   `DAE_FORTRAN_STYLE` - If ON, the matrices will be defined using FORTRAN style (one-based indexing of columns and rows). By default it is OFF (zero-based indexing).
-   `DAE_SINGLE` - If ON, the single precision will be used in the solver instead of double. Single precision may ruin the accuracy. It is highly recommended to leave this option OFF. This option exists for the future compatibility with CUDA implementations of the solver.
-   `DAE_BUILD_EXAMPLES` - Build all the examples, ON by default.
-   `DAE_TEST` - Build automatic solver test, ON by default. The test can be executed by the command `ctest` from the building directory.
-   `DAE_MKL_DIR` - Defines a path to Intel MKL root directory (usually `/opt/intel/mkl`).
-   `PLOTTING` - Use [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module for plotting, `OFF` by default. If `ON`, `cmake` will try to find Python and `numpy` include directories and libraries.
-   `PYTHON_INCLUDE` - Only if `PLOTTING=ON`, defines a path to Python include file (`Python.h`) for plotting.
-   `PYTHON_NUMPY_INCLUDE` - Only if `PLOTTING=ON`, defines a path to Python `numpy` include file (`numpy/arrayobject.h`) for plotting.
-   `PYTHON_LIB` - Only if `PLOTTING=ON`, defines Python library (e.g. `libpython3.6m`) for plotting.

### Windows

Download and install compiler (e.g. [Microsoft Visual Studio](https://visualstudio.microsoft.com/downloads/)) and [Python 3](https://www.python.org/downloads/windows/) with `numpy` and `matplotlib` modules (for plotting, optional).

Download and install [Git](https://git-scm.com/download/win) and [CMake](https://cmake.org/download/) for Windows.

From `Git Bash` command line clone dae-cpp library (you may need to create a working directory first):

```bash
git clone https://github.com/ikorotkin/dae-cpp.git
```

Start CMake (`cmake-gui`), choose the source code path (`dae-cpp` folder) and empty target directory (it will contain Visual Studio project files). Press "Configure" button.

If CMake cannot find any of the libraries, it will print an error message. You can modify the paths and other parameters (see [More building options](https://github.com/ikorotkin/dae-cpp#more-building-options) above) and re-configure the project.

If configuration is successful, press "Configure" again to update the cache and then "Generate". In the target directory you will find Visual Studio project files.

Double-click on `dae-cpp.sln` to open Visual Studio with the project. Do not forget to change Solution Configuration from `Debug` to `Release`. Build the solution (`F7` by default). After compilation, the executable files can be found in `Release` folder.

Note that in order to execute the tests (for example, `perovskite.exe`) from `Release` folder, you need to set up Intel MKL environment variables by executing `mklvars.bat intel64` or `mklvars.bat ia32` (depending on the target platform) from `cmd`. By default `mklvars.bat` is located in MKL root folder in `bin` subdirectory, for example:

```bash
"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64
```

**_Alternatively_**, you may install [Windows Subsystem for Linux](https://docs.microsoft.com/en-gb/windows/wsl/install-win10?redirectedfrom=MSDN) and your preferred Linux Distribution (e.g. Ubuntu), and then just follow [installation instructions for Linux](#linux).

### Mac

Make sure you have installed: `git`, `cmake`, `gcc` and, optional, Python 3 with `numpy` and `matplotlib` modules (for plotting). If these packages are not installed yet, you may install [Homebrew](https://brew.sh/) (package manager for macOS), then install all necessary packages:

```bash
brew install cmake git gcc python
pip3 install numpy matplotlib
```

Note if you install `git` for the first time you will need to configure it (change `Your Name` and `your@email` to your full name and email):

```bash
git config --global user.name "Your Name"
git config --global user.email your@email
```

Then from the working directory download dae-cpp library source files:

```bash
git clone https://github.com/ikorotkin/dae-cpp.git
```

Check the version of `gcc` compiler by typing `gcc` and pressing `Tab` key a few times in the terminal, it will show you the version of `gcc` currently installed, for example, `gcc-9` (you could use the command `gcc --version` but it may point to `clang` compiler for Mac that does not support OpenMP out of the box).

Create `build` directory:

```bash
cd dae-cpp
mkdir build
cd build
```

Configure the project. *Make sure `g++` version (9 in the example below) is correct*:

```bash
cmake .. -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX=$PWD
```

In the command above you may change the user-defined path where the package should be installed (type it instead of `$PWD`). By default the package will be installed into the current `build` directory.

Note that `cmake` will try to find Intel MKL at its default location `/opt/intel/mkl` or according to `MKLROOT` environment variable. If the installation path is different, please provide MKL root path with the following `cmake` option: `-DDAE_MKL_DIR=/path_to_intel_mkl`.

Instead of `cmake ..` you may consider using `ccmake ..`, a UI for `cmake` that will allow you to see and change all the options available before building the solver.

Install dae-cpp and perform a quick self test:

```bash
make -j4
make install
ctest
```

## How to use

Please refer to [simple_dae.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/simple_dae/simple_dae.cpp), [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp), [perovskite.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.cpp) or [diffusion_2d.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d.cpp) as an example.

The main usage algorithm can be the following. Consider we have a system of DAEs written in a matrix-vector form, with some Mass matrix, RHS, and some initial conditions.

### Step 0. Include dae-cpp into the project

Include the solver's main header to the project. A shortcut to the solver's namespace (`daecpp`) can be added as well:

```cpp
#include "path/to/dae-cpp/include/solver.h"
namespace dae = daecpp;
```

### Step 1. Define the DAE parameters and initial state vector

For example, for *N* equations we should define the state vector with the size *N* and initialize it in accordance with the initial conditions:

```cpp
// State vector
dae::state_type x(N);

// Initial conditions
for(MKL_INT i = 0; i < N; i++)
{
    x[i] = 1.0;
}
```

We can get access to each element of the state vector **x** as to `std::vector` from STL. Also note that instead of `int` or any other integer types we should use `MKL_INT` type. This gives us possibility to re-compile the project with `DAE_LONG_INT` option, so the code will work fine even for extremely huge systems (with *N* more than 10<sup>7</sup>).

### Step 2. Set up the RHS

Create MyRHS class that inherits the abstract `daecpp::RHS` class from dae-cpp library. The parent RHS class contains a pure virtual functor (operator `()`), that must be overridden in the child class. See, for example, [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp), [perovskite_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_RHS.cpp) or [diffusion_2d_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d_RHS.cpp).

Once the RHS class is overridden, we can create an instance of the child class with some user-defined parameter container *p*:

```cpp
MyRHS rhs(p);
```

In the child MyRHS class the user can also override `stop_condition` virtual function. By default (if not overridden) the function always returns `false`. The user may override this behaviour and set up one or several stop conditions for the solver depending on the solution **x** at the current time *t*. As soon as the function returns `true`, the solver will finalise the current time step and return the current solution. A trivial example of the stop condition function can be found in [perovskite_RHS.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_RHS.cpp).

For the debugging purposes, the RHS can be saved to a file:

```cpp
rhs.dump(x, 0);
rhs.dump(x, 0.1);
```

In this example we saved two RHS vectors, at time 0 and 0.1.

### Step 3. Set up the Mass matrix

Create MyMassMatrix class that inherits the abstract `daecpp::MassMatrix` class from dae-cpp library. Similar to the previous step, the parent MassMatrix class contains a pure virtual functor (operator `()`), that must be overridden in the child class. Refer to [perovskite_Mass.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Mass.cpp) or [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp) as an example. Note that the matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format). See also [a note about Sparse Matrix Format](#a-note-about-Sparse-Matrix-Format).

Create an instance of the child MyMassMatrix class with the given size *N*:

```cpp
MyMassMatrix mass(N);
```

If the Mass matrix is a simple identity matrix, one can use `daecpp::MassMatrixIdentity` class from dae-cpp library instead of inheriting `daecpp::MassMatrix`. This will create identity Mass matrix in sparse format with the given size *N*:

```cpp
dae::MassMatrixIdentity mass(N);
```

For the debugging purposes, you can save the Mass matrix to a file:

```cpp
mass.dump();
```

### Step 4. Set up Jacobian matrix

We can provide analytical Jacobian by overriding `daecpp::Jacobian` class from the dae-cpp library (see [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp) or [perovskite_Jacobian.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_Jacobian.cpp)) or just use numerically estimated one (this may significantly slow down the computation for large *N*). If provided, analytical Jacobian matrix should be defined in [three array sparse format](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format) similar to the Mass matrix. See also [a note about Sparse Matrix Format](#a-note-about-Sparse-Matrix-Format).

If we don't provide analytical Jacobian we should estimate it with the given tolerance:

```cpp
dae::Jacobian jac(rhs, 1.0e-6);
```

Note that we should pass an instance of the user-defined RHS in order to estimate numerical Jacobian.

Again, for the debugging purposes, Jacobian can be saved to a file:

```cpp
jac.dump(x, 0);
jac.dump(x, 0.1);
```

In the example above we saved two Jacobians, at time 0 and 0.1.

In some cases the derivation and coding of the analytic Jacobian can be a tricky problem itself. So `dae::Jacobian` class provides additional functionality to compare two Jacobians (one of them is numerical) and write the differences:

```cpp
dae::Jacobian jac(rhs, 1.0e-6);  // Numerical Jacobian calculated automatically (slow)
MyJacobian jac_user(rhs);        // Analytic Jacobian provided by the user

// Comparison of jac and jac_user and writing the differences to a file
jac_user.compare(jac, x, 0.1, 1e-4);
```

Here we compared two Jacobians at time 0.1 with the relative tolerance 10<sup>-4</sup>.

### Step 5. Set the solver options

The solver has lots of options related to the solution process. They all have some default values (defined in [solver_options.h](https://github.com/ikorotkin/dae-cpp/blob/master/src/solver_options.h)) but they can be overridden by a user:

```cpp
// Create an instance of the solver options and update some of the solver
// parameters defined in solver_options.h
dae::SolverOptions opt;

// For example, let's change the initial time step
opt.dt_init = 0.01;
```

### Step 6. Solve the system

Now we are ready to create an instance of the solver with particular RHS, Mass matrix, Jacobian and the solver options, and then start the solver:

```cpp
dae::Solver solve(rhs, jac, mass, opt);
int status = solve(x, t1);
```

Here *t*<sub>1</sub> is the integration time (0 < *t* < *t*<sub>1</sub>), and **x** is the initial condition vector defined above.

The solver returns 0 if integration is successful or error code otherwise. Solution at time *t*<sub>1</sub> will be written into vector **x** (initial conditions will be overwritten). The actual integration time *t*<sub>1</sub> will be returned (the solver may terminate integration earlier). That's it!

#### Optional: Set up Observer

In order to get intermediate solutions at times *t*<sub>a</sub>, *t*<sub>b</sub>, *t*<sub>c</sub>, etc. (0 < *t*<sub>a</sub> < *t*<sub>b</sub> < *t*<sub>c</sub> < ... < *t*<sub>1</sub>), for example, for plotting, one can call the solver at the given times:

```cpp
solve(x, t_a);  // solves the system in the interval [0; t_a] and stores the solution in x
solve(x, t_b);  // continues solving in the interval [t_a; t_b], replaces the solution in x
solve(x, t_c);  // continues solving in the interval [t_b; t_c], replaces the solution in x
                // ...
solve(x, t1);   // continues solving in the interval [t_c; t1] and
                // stores the final solution at time t1 in x
```

Every call the solver will take the previous solution **x** (if available from the previous call) and overwrite it with a new one at the given time.

But a proper (and more efficient) way to get intermediate results is to override `virtual void observer(...)` function from `daecpp::Solver` class. This observer function receives the current solution vector **x** and the current time *t* every time step and allows a user to get access to the solution at each time layer. An example of a simple observer is given in the file [robertson.cpp](https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.cpp), also in [perovskite_observer.h](https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite_observer.h).

### Step 7 (optional). Plot results

Solution can be visualised using a simple [C++ interface](https://github.com/lava/matplotlib-cpp) to Python [matplotlib](https://matplotlib.org/) module. For example, if `python`, `numpy` and `matplotlib` are installed and the solver was built with `PLOTTING=ON`, the [perovskite](https://github.com/ikorotkin/dae-cpp/tree/master/examples/perovskite) example will produce the following plot:

<p align="center">
  <img src="https://github.com/ikorotkin/dae-cpp/blob/master/examples/perovskite/perovskite.png">
</p>

Here *P(x)* is the ion concentration in a perovskite solar cell, and *Phi(x)* is the corresponding potential distribution.

The second example, [diffusion_2d](https://github.com/ikorotkin/dae-cpp/tree/master/examples/diffusion_2d), will produce a two-dimensional Gaussian function, a solution of two-dimensional diffusion problem with an instantaneous point source in the middle of the plane:

<p align="center">
  <img src="https://github.com/ikorotkin/dae-cpp/blob/master/examples/diffusion_2d/diffusion_2d.png">
</p>

The third example, [robertson](https://github.com/ikorotkin/dae-cpp/tree/master/examples/robertson), solves [Robertson stiff DAE problem](https://www.mathworks.com/help/matlab/ref/ode15s.html) with a conservation law. It produces the following figure:

<p align="center">
  <img src="https://github.com/ikorotkin/dae-cpp/blob/master/examples/robertson/robertson.png">
</p>

Note that by default the plotting is switched off in the examples, but the plotting-related code can be activated using `#define PLOTTING` at the very beginning of each example. Activating the plotting refers to `matplotlibcpp.h` header located in `src/external/matplotlib-cpp/` directory.

### A note about Sparse Matrix Format

It should be noted that you must define all the diagonal elements of the matrix, even if they are zero. This greatly increases performance, and if some rows are skipped, the code will just stop working. Please double check your Mass matrix and Jacobian, they both should have the main diagonal filled in. Even if the given row is empty (all elements are zero), define zero on the main diagonal explicitly.

If you are struggling with Intel MKL sparse format, you can use simple three-array format instead, where you need to define all non-zero elements and their indexes (coordinates) in the matrix. For example for the identity 3x3 matrix, you only need to define three non-zero elements and their position in the matrix:

```cpp
M.A.resize(3);   // Number of non-zero elements
M.ia.resize(3);  // Number of non-zero elements
M.ja.resize(3);  // Number of non-zero elements

M.A[0] = 1;   // First non-zero or diagonal element
M.ia[0] = 0;  // Column index of the first non-zero element
M.ja[0] = 0;  // Raw index of the first non-zero element

M.A[1] = 1;   // Second non-zero or diagonal element
M.ia[1] = 1;  // Column index of the second non-zero element
M.ja[1] = 1;  // Raw index of the second non-zero element

M.A[2] = 1;   // Third non-zero or diagonal element
M.ia[2] = 2;  // Column index of the third non-zero element
M.ja[2] = 2;  // Raw index of the third non-zero element
```

This form will be automatically converted to the three-array sparse format compatible with Intel MKL. Do not forget to define all diagonal elements even if they are zero. Do not mix the elements up (fill in the first row from left to right, then the second row, etc.).

## Contribution and feedback

Please feel free to contribute into the project!

If you have any questions, suggestion, or a feedback, please, submit an [issue](https://github.com/ikorotkin/dae-cpp/issues).

## Licensing

-   dae-cpp is fully open source under [MIT license](https://github.com/ikorotkin/dae-cpp/blob/master/LICENSE).
-   Intel MKL is free for use and redistribution under [Intel Simplified Software License](https://software.intel.com/en-us/license/intel-simplified-software-license).
