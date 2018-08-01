#include <stdio.h>

#include "molecule.h"


int main(){
    FILE *molxyz;
    Molecule mol;

    molxyz =  fopen("acetaldehyde.txt", "r");
    molecule_read(molxyz, &mol);

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
    return 0;
}

