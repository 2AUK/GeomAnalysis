#include <stdio.h>

#include "molecule.h"


int main(){
    FILE *molxyz;
    Molecule mol;

    molxyz =  fopen("acetaldehyde.txt", "r");
    molecule_read(molxyz, &mol);


    printf("%d\n", mol.natom);
    
    return 0;
}

