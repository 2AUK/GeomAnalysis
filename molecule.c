#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"
#include "LAS.h"
#include "masses.h"
    
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

    //Need to allocate a array for mol.geom in contiguous memory
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

double molecule_oop(Molecule mol, int a, int b, int c, int d){
    double ebcd_x = (unit_vector(mol, 1, c, b) * unit_vector(mol, 2, c, d) - unit_vector(mol, 2, c, b) * unit_vector(mol, 1, c, d));
    double ebcd_y = (unit_vector(mol, 2, c, b) * unit_vector(mol, 0, c, d) - unit_vector(mol, 0, c, b) * unit_vector(mol, 2, c, d));
    double ebcd_z =  (unit_vector(mol, 0, c, b) * unit_vector(mol, 1, c, d) - unit_vector(mol, 1, c, b) * unit_vector(mol, 0, c, d));

    double exx = ebcd_x * unit_vector(mol, 0, c, a);
    double eyy = ebcd_y * unit_vector(mol, 1, c, a);
    double ezz = ebcd_z * unit_vector(mol, 2, c, a);

    double theta = (exx + eyy + ezz) / sin(molecule_angle(mol, b, c, d));

    if (theta < -1.0) theta = asin(-1.0);
    else if (theta > 1.0) theta = asin(1.0);
    else theta =  asin(theta);

    return theta;
}

double molecule_torsion(Molecule mol, int a, int b, int c, int d){
    double eabc_x = (unit_vector(mol,1,b,a) * unit_vector(mol,2,b,c) - unit_vector(mol,2,b,a)*unit_vector(mol,1,b,c));
    double eabc_y = (unit_vector(mol,2,b,a) * unit_vector(mol,0,b,c) - unit_vector(mol,0,b,a)*unit_vector(mol,2,b,c));
    double eabc_z = (unit_vector(mol,0,b,a) * unit_vector(mol,1,b,c) - unit_vector(mol,1,b,a)*unit_vector(mol,0,b,c));

    double ebcd_x = (unit_vector(mol,1,c,b) * unit_vector(mol,2,c,d) - unit_vector(mol,2,c,b)*unit_vector(mol,1,c,d));
    double ebcd_y = (unit_vector(mol,2,c,b) * unit_vector(mol,0,c,d) - unit_vector(mol,0,c,b)*unit_vector(mol,2,c,d));
    double ebcd_z = (unit_vector(mol,0,c,b) * unit_vector(mol,1,c,d) - unit_vector(mol,1,c,b)*unit_vector(mol,0,c,d));

    double exx = eabc_x * ebcd_x;
    double eyy = eabc_y * ebcd_y;
    double ezz = eabc_z * ebcd_z;

    double tau = (exx + eyy + ezz) / (sin(molecule_angle(mol, a, b, c)) *  sin(molecule_angle(mol, b, c, d)));

    if (tau < -1.0) tau = acos(-1.0);
    else if (tau > 1.0) tau = acos(1.0);
    else tau = acos(tau);

    double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
    double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
    double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
    double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;
    double sign = 1.0;
    double dot = cross_x*unit_vector(mol, 0, b, c)+cross_y*unit_vector(mol, 1, b, c)+cross_z*unit_vector(mol,2,b,c);
    if (dot < 0.0) sign = -1.0;

    return tau*sign;
}

void molecule_COMtranslation(Molecule *mol){
    double M = 0.0;
    
    for (int i = 0; i < mol->natom; i++) M += masses[mol->zvals[i]];

    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    double mi;

    for (int i = 0; i < 3*mol->natom; i+=3){
        mi = masses[mol->zvals[i/3]];
        xcm += mi * mol->geom[i+0];
        ycm += mi * mol->geom[i+1];
        zcm += mi * mol->geom[i+2];
    }

    xcm /= M;
    ycm /= M;
    zcm /= M;

    printf("COM is\t%lf\t%lf\t%lf\n", xcm, ycm, zcm);

    for (int i = 0; i < 3*mol->natom; i+=3){
        mol->geom[i+0] += -xcm;
        mol->geom[i+1] += -ycm;
        mol->geom[i+2] += -zcm;
    }
}



void molecule_free(Molecule *mol){
    free(mol->zvals);
    free(mol->geom);
}
