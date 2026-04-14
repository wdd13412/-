# AGENTS.md

## Cursor Cloud specific instructions

### Project overview

This is a Fortran-based 2D compressible Euler CFD solver (BuFlow) with Tapenade-generated automatic differentiation for adjoint-based aerodynamic shape optimization of an RAE 2822 airfoil. There are no web services, databases, or package managers — just Fortran source files compiled into a single executable.

### System dependencies

- **gfortran** (GNU Fortran compiler, v13+)
- **LAPACK/BLAS** dev libraries (`liblapack-dev`, `libblas-dev`)

Install with: `sudo apt-get install -y gfortran liblapack-dev libblas-dev`

### Build

Compile from the workspace root (order matters for module dependencies):

```bash
gfortran -O2 -o cfd_solver TypesModule_d.f90 meshdeformationn_d.f90 BuFlow_test_d.f90 main_d.f90 run_parameter_d.f90 -llapack -lblas
```

### Run

The solver must be executed from `/workspace` since it reads mesh data and CSV files via relative paths:

```bash
cd /workspace && ./cfd_solver
```

Execution takes approximately 3–4 minutes. Output:
- `solution_001.vtk` — VTK visualization file (viewable in ParaView)
- `wall_cp_data.csv` — wall pressure coefficient distribution with tangent-mode gradients
- `wing_update.csv` — deformed airfoil coordinates

### Key caveats

- **No build system**: There is no Makefile or CMakeLists.txt. The compilation command above must be used directly.
- **No tests**: There is no automated test suite. Verification is done by running the solver and checking that it produces output files and prints "CFD计算完成" (CFD computation complete) at the end.
- **No lint**: There are no Fortran linters configured. Static analysis can be done with `gfortran -Wall -Wextra -fsyntax-only *.f90` for syntax checking.
- **Module files**: Compilation generates `.mod` files in the working directory. These are build artifacts and should not be committed.
- **Long runtime**: The solver runs ~18,500 cells through iterative time stepping and GMRES linear solver, taking 3–4 minutes to complete.
