## Reproducible build for `true||r||/||b|| ~= 0.1324`

Use the exact command below (note the `-O2`):

`gfortran -O2 -o buflow_test_d_opt TypesModule_d.f90 meshdeformationn_d.f90 BuFlow_test_d.f90 main_d.f90 run_parameter_d.f90 -llapack -lblas`

Then run:

`./buflow_test_d_opt`

Expected final GMRES line at `outer=10`:

`true||r||/||b||=  0.13240372991858221`

## Important

- If you compile **without** `-O2`, the same source gives about `0.16509053348880676`.
- `BuFlow_test_d.f90` and `BuFlow_test_d_final.f90` are byte-identical in this branch.

SHA256 (both files):

`012d41ec95362b2a1372c350b96f2e82fe7383ff8fee3197825ab2c1c11815b1`
