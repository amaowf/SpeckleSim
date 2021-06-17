#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

mkdir --parent src/lib/tree
mkdir --parent src/lib/mlfmm/config


ifort -O0 -g -c -openmp -cpp -o "src/lib/lib_hash_function.o" "../src/lib/lib_hash_function.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/lib_sort.o" "../src/lib/lib_sort.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/toolbox.o" "../src/lib/toolbox.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/file_io.o" "../src/lib/file_io.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/tree/lib_tree_type.o" "../src/lib/tree/lib_tree_type.f90"
ifort -O0 -g -c -cpp -v -o "src/lib/tree/lib_tree_helper_functions.o" "../src/lib/tree/lib_tree_helper_functions.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/tree/lib_tree_type_operator.o" "../src/lib/tree/lib_tree_type_operator.f90"
ifort -O0 -g -c -cpp -o "src/lib/tree/lib_tree.o" "../src/lib/tree/lib_tree.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/tree/lib_tree_public.o" "../src/lib/tree/lib_tree_public.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/mlfmm/config/ml_fmm_type.o" "../src/lib/mlfmm/config/ml_fmm_type.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/mlfmm/lib_ml_fmm_type_operator.o" "../src/lib/mlfmm/lib_ml_fmm_type_operator.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/mlfmm/config/ml_fmm_math.o" "../src/lib/mlfmm/config/ml_fmm_math.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/mlfmm/lib_ml_fmm_helper_functions.o" "../src/lib/mlfmm/lib_ml_fmm_helper_functions.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/mlfmm/lib_ml_fmm.o" "../src/lib/mlfmm/lib_ml_fmm.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/light _scattering_by_particles/light_scattering.o" "../src/lib/light _scattering_by_particles/light_scattering.f90"
ifort -O0 -g -c -openmp -cpp -o "src/lib/lib_data_types.o" "../src/lib/lib_data_types.f90"
ifort -O0 -g -c -openmp -cpp -o "src/main.o" "../src/main.f90"
ifort -openmp -o "FoggySim"  ./src/lib/tree/lib_tree.o ./src/lib/tree/lib_tree_helper_functions.o ./src/lib/tree/lib_tree_public.o ./src/lib/tree/lib_tree_type.o ./src/lib/tree/lib_tree_type_operator.o  ./src/lib/mlfmm/config/ml_fmm_math.o ./src/lib/mlfmm/config/ml_fmm_type.o  ./src/lib/mlfmm/lib_ml_fmm.o ./src/lib/mlfmm/lib_ml_fmm_helper_functions.o ./src/lib/mlfmm/lib_ml_fmm_type_operator.o  ./src/lib/light\ _scattering_by_particles/light_scattering.o  ./src/lib/file_io.o ./src/lib/lib_data_types.o ./src/lib/lib_hash_function.o ./src/lib/lib_sort.o ./src/lib/toolbox.o  ./src/main.o  ./.metadata/.plugins/org.eclipse.cdt.make.core/specs.o   