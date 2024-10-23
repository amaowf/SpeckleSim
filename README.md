SpeckleSimulator was developed at ITO, University of Stuttgart uisng boundary element method (or surface integral equation method). Two meshing methods are implemented, one is flat triangular element and the other one is quadrilateral higher order 10 edges element. 

To treat the lower and higher order sigularity problem, the methods in this paper were implemented: I. Hänninen, M. Taskinen, and J. Sarvas, "Singularity subtraction integral formulae for surface integral equations with RWG, rooftop and hybrid basis functions," Progress in Electromagnetics Research 63, 243–278 (2006).

GMRES was used to solve the linear equation and the fortran code was taken from: A set of Flexible-GMRES routines for real and complex arithmetics (https://www.cerfacs.fr/algor/reports/1998/TR_PA_98_20.pdf)

