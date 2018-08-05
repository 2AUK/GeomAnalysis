#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

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

void molecule_inertia(Molecule mol){
    gsl_matrix * inertia_mat = gsl_matrix_alloc (3, 3);
    double mi, I00, I11, I22, I01, I02, I12;
    I00 = I11 = I22 = I01 = I02 = I12 = 0;

    for (int i = 0; i < 3*mol.natom; i+=3){
        mi = masses[mol.zvals[i/3]];
        I00 += mi * (mol.geom[i+1]*mol.geom[i+1] + mol.geom[i+2]*mol.geom[i+2]);
        gsl_matrix_set(inertia_mat, 0, 0, I00);
        I11 += mi * (mol.geom[i+0]*mol.geom[i+0] + mol.geom[i+2]*mol.geom[i+2]);
        gsl_matrix_set(inertia_mat, 1, 1, I11);
        I22 += mi * (mol.geom[i+0]*mol.geom[i+0] + mol.geom[i+1]*mol.geom[i+1]);
        gsl_matrix_set(inertia_mat, 2, 2, I22);
        I01 += mi * (mol.geom[i+0]*mol.geom[i+1]); 
        gsl_matrix_set(inertia_mat, 0, 1, I01);
        I02 += mi * (mol.geom[i+0]*mol.geom[i+2]); 
        gsl_matrix_set(inertia_mat, 0, 2, I02);
        I12 += mi * (mol.geom[i+1]*mol.geom[i+2]); 
        gsl_matrix_set(inertia_mat, 1, 2, I12);
    }

    gsl_matrix_set(inertia_mat, 1, 0, gsl_matrix_get(inertia_mat, 0, 1));
    gsl_matrix_set(inertia_mat, 2, 0, gsl_matrix_get(inertia_mat, 0, 2));
    gsl_matrix_set(inertia_mat, 2, 1, gsl_matrix_get(inertia_mat, 1, 2));

    printf("++++++Inertia Tensor (amu bohr^2)++++++\n");
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++)
            printf("%lf\t", gsl_matrix_get(inertia_mat, i, j));
        putchar('\n');
    }

    gsl_eigen_symmv_workspace *eig_solver = gsl_eigen_symmv_alloc(3);
    gsl_vector *evals = gsl_vector_alloc(3);
    gsl_matrix *evecs = gsl_matrix_alloc(3, 3);

    gsl_eigen_symmv(inertia_mat, evals, evecs, eig_solver);
    
    
    printf("++++++Principal Moments of Inertia (amu bohr^2)++++++\n");
    for(int i = 0; i< 3; i++)
        printf("%lf\n", gsl_vector_get(evals, i));
    double conv1 = 0.529177249 * 0.529177249;
    printf("++++++Principal Moments of Inertia (amu  AA^2)++++++\n");
    for(int i = 0; i<3; i++)
        printf("%lf\n", gsl_vector_get(evals, i) * conv1);
    double conv2 = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
    printf("++++++Principal Moments of Inertia (g cm^2)++++++\n");
    for(int i = 0; i<3; i++)
        printf("%lf\n", gsl_vector_get(evals, i) * conv2);
    if (mol.natom == 2)  printf("Molecule is diatomic\n");
    else if(gsl_vector_get(evals, 0) < 1E-4) printf("Molecule is linear\n");
    else if(fabs(gsl_vector_get(evals, 0) - gsl_vector_get(evals, 1)) < 1E-4 && fabs(gsl_vector_get(evals, 1) - gsl_vector_get(evals, 2)) < 1E-4)
        printf("Molecule is a spherical top\n");
    else if(fabs(gsl_vector_get(evals, 0) - gsl_vector_get(evals, 1)) < 1E-4 && fabs(gsl_vector_get(evals, 1) - gsl_vector_get(evals, 2)) > 1E-4)
        printf("Molecule is an oblate spherical top\n");
    else if(fabs(gsl_vector_get(evals, 0) - gsl_vector_get(evals, 1)) > 1E-4 && fabs(gsl_vector_get(evals, 1) - gsl_vector_get(evals, 2)) < 1E-4)
        printf("Molecule is a prolate spherical top\n");
    else printf("Molecule is an asymmetric top\n");


    double _pi = acos(-1.0);
    double conv = 6.6260755E-34/(8.0 * _pi * _pi);
    conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
    conv *= 1E-6;
    printf("++++++Rotational Constants (MHz)++++++\n");
    printf("\tA = %lf\tB = %lf\tC = %lf\n", conv/gsl_vector_get(evals, 0), conv/gsl_vector_get(evals, 1), conv/gsl_vector_get(evals, 2));

    double convc = 6.6260755E-34 / (8.0 * _pi * _pi);
    convc /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
    convc /= 2.99792458E10;
    printf("++++++Rotational Constants (cm^-1)++++++\n");
    printf("\tA = %lf\tB = %lf\tC = %lf\n", convc/gsl_vector_get(evals, 0), convc/gsl_vector_get(evals, 1), convc/gsl_vector_get(evals, 2));


    gsl_vector_free(evals);
    gsl_matrix_free(evecs);
    gsl_matrix_free(inertia_mat);
    gsl_eigen_symmv_free(eig_solver);
}
    
void molecule_free(Molecule *mol){
    free(mol->zvals);
    free(mol->geom);
}
