#include <stdio.h>

#include "molecule.h"


int main(){
    FILE *molxyz;
    Molecule mol;

    molxyz =  fopen("acetaldehyde.txt", "r");
    molecule_read(molxyz, &mol);
    printf("++++++Input Data++++++");
    printf("%d\n", mol.natom);
    for (int i = 0; i < mol.natom; i++)
        printf("%d\t", mol.zvals[i]);
    putchar('\n');
    for (int j = 0; j < 3*mol.natom; j+=3){
        for (int k = 0; k < 3; k++){
            printf("%3.6lf\t", mol.geom[j+k]);
        }
        putchar('\n');
    }
    printf("+++++++Distances between Atoms+++++++\n");
    for (int n = 0; n < 3*mol.natom; n+=3){
        for (int m = 0; m < n; m+=3){
            printf("%d\t%d\t%lf\n", n/3, m/3, molecule_bond(mol, n, m));
        }
    }
    return 0;
}

