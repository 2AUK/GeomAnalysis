#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"
#include "LAS.h"
    
#define BUF_SIZE 20000

void molecule_read(FILE *stream, Molecule *mol){
    //Format of .xyz file
    // 7                                      - Number of atoms in molecule
    // 16 0.000 0.000 0.000                   - z-value x y z
    // 1        ...
    // ...
    // ...
    //And so on for however many atoms are in the molecule
    
    char c;
    char buf[BUF_SIZE];
    char input_natom[3];
    int n, i, j, count;
    n = i = j = count = 0;
    double x, y, z, z_val;

    //Reading number of atoms   
    while((c = fgetc(stream)) != '\n')
        input_natom[count++] = c;
    mol->natom = atoi(input_natom);

    //Need to allocate a 2x2 by array for mol.geom in contiguous memory
    mol->zvals = malloc(sizeof(int) * mol->natom);
    mol->geom = malloc(sizeof(double) * mol->natom * mol->natom);

    //Read and allocate z-values and coordinates in to respective arrays
    while (fgets(buf, BUF_SIZE, stream) != NULL){
        sscanf(buf, "%lf %lf %lf %lf", &z_val, &x, &y, &z);
        mol->zvals[i++] = z_val;
        mol->geom[j++] = x;
        mol->geom[j++] = y;
        mol->geom[j++] = z;
    }
}

double molecule_bond(Molecule mol, int a, int b){
    return sqrt( (mol.geom[a+0] - mol.geom[b+0]) * (mol.geom[a+0] - mol.geom[b+0])
               + (mol.geom[a+1] - mol.geom[b+1]) * (mol.geom[a+1] - mol.geom[b+1])
               + (mol.geom[a+2] - mol.geom[b+2]) * (mol.geom[a+2] - mol.geom[b+2]));
}

double unit_vector(Molecule mol, int cart, int a, int b){
     return -1 * ((mol.geom[a + cart] - mol.geom[b + cart]) / molecule_bond(mol, a, b));    
}
double molecule_angle(Molecule mol, int a, int b, int c){
    return acos(unit_vector(mol, 0, b, a) * unit_vector(mol, 0, b, c) + unit_vector(mol ,1, b, a) * unit_vector(mol, 1, b, c) + unit_vector(mol, 2, b, a) * unit_vector(mol, 2, b, c));
}


void molecule_free(Molecule *mol){
    free(mol->zvals);
    free(mol->geom);
}
