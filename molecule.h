#ifndef HEADERFILE_H
#define HEADERFILE_H
#include <stdio.h>

// Molecule struct to hold all information on molecule
typedef struct Molecule{
    int natom;
    int charge;
    int *zvals;
    double *geom;
}Molecule;

// Functions related to initialisation and manipulation and freeing of Molecule struct

void molecule_read(FILE *stream, Molecule *mol);

double molecule_bond(Molecule mol, int a, int b);

double unit_vector(Molecule mol, int cart, int a, int b);

double molecule_angle(Molecule mol, int a, int b, int c);

double molecule_oop(Molecule mol, int a, int b, int c, int d);

double molecule_torsion(Molecule mol, int a, int b, int c, int d);

double molecule_COMtranslation(Molecule *mol);

void molecule_free(Molecule *mol);
#endif
    
