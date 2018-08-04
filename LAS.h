#ifndef HEADERFILE_H
#define HEADERFILE_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "molecule.h"

//helper functions for performing intermediate mathematical operations
double unit_vector(Molecule mol, int cart, int a, int b);

//wrappers around the annoyingly complicated gsl routines TODO
#endif
