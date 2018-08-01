#ifndef HEADERFILE_H
#define HEADERFILE_H
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

#include "molecule.h"

//helper functions for performing intermediate mathematical operations
void unit_vector(Molecule mol_i, Molecule mol_j);

//wrappers around the annoyingly complicated gsl routines TODO
#endif
