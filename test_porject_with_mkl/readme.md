# How to use the Intel MKL

## Install the Intel MKL

- Download the [Intel MKL](https://software.intel.com/en-us/mkl).
- Install the Intel MKL

## Setup of the environment variable *MKLROOT*

Operating System:
- Linux (Ubuntu)
 
```bash
$ echo export MKLROOT=<path-to>/intel/mkl >> ~/.bashrc
```
   
- Windows
    - Open the environment variable editor and add the *MKLROOT* variable
    - Set the value of the variable to the path of the intel mkl folder
    
       MKLROOT: <path-to>/intel/mkl
   

## Hint
If other libraries from the Intel MKL are required, use the Intel Math Kernel Library Link Line Advisor to customize the project properties.

[Webpage: Link Line Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) 