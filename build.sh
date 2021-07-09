#!/bin/bash

cpu_cores=$(nproc --all)
FoBiS.py build -fc "gfortran-8 -cpp -fopenmp -ffree-line-length-none -frecursive" -j $cpu_cores -ed FoggySim src/lib/libmath_generated_functions src/lib/libmath/src/math/lib_math_wigner.f90 .ipynb_checkpoints
