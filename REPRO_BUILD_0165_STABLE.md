## Why you get `0.165` with same source

This project is numerically sensitive to compiler floating-point transformations.
With the same `BuFlow_test_d.f90`, different optimization/vectorization choices
lead to different GMRES trajectories.

Observed in this branch:

- `gfortran -o buflow_test_d ... -llapack` -> `true||r||/||b|| ~= 0.1650905335`
- `gfortran -O2 -fno-tree-vectorize -o buflow_test_d ... -llapack -lblas` -> `~0.1650905335`
- `gfortran -O2 -o buflow_test_d ... -llapack -lblas` -> `~0.1324037299`

So this is not caused by missing packages or skipped code paths; it is floating-point
execution-order sensitivity.

## Stable command to match your Ubuntu behavior

Use this command if you want robustly reproducible `~0.165` across machines:

`gfortran -O2 -fno-tree-vectorize -o buflow_test_d TypesModule_d.f90 meshdeformationn_d.f90 BuFlow_test_d.f90 main_d.f90 run_parameter_d.f90 -llapack -lblas`

Then run:

`./buflow_test_d`
