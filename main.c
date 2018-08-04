#include <stdio.h>
#include <math.h>

#include "molecule.h"

int main(){
    FILE *molxyz;
    Molecule mol;

    molxyz =  fopen("acetaldehyde.txt", "r");
    molecule_read(molxyz, &mol);
    printf("++++++Input Data++++++\n");
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
    printf("+++++++Bond Angles+++++++\n");
    for (int a = 0; a < 3*mol.natom; a+=3){
        for (int s = 0; s < a; s+=3){
            for (int d = 0; d < s; d+=3){
                if (molecule_bond(mol, a, s) < 4.0 && molecule_bond(mol, s, d) < 4.0)
                    printf("%d-%d-%d\t%lf\n", a/3, s/3, d/3, molecule_angle(mol, a, s, d)*(180.0/acos(-1.0)));
            }
        }
    }
    return 0;
}

