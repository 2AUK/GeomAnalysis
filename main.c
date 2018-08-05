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
    printf("++++++Out Of Plane Angles++++++\n");
    for (int i = 0; i < 3*mol.natom; i+=3){
        for (int k = 0; k < 3*mol.natom; k+=3){
            for (int j = 0; j < 3*mol.natom; j+=3){
                for(int l = 0; l < j; l+=3){
                    for( int count = 0; count < 3; count++){
                        if (i+count != j+count
                          &&i+count != k+count
                          &&i+count != l+count
                          &&j+count != k+count
                          &&j+count != l+count
                          &&k+count != l+count
                          &&molecule_bond(mol, i, k) < 4.0
                          &&molecule_bond(mol, k, j) < 4.0
                          &&molecule_bond(mol, k, l) < 4.0){
                            if (count == 0)
                                printf("%d-%d-%d-%d\t%lf\n", i/3, j/3, k/3, l/3, molecule_oop(mol,i,j,k,l)*(180.0/acos(-1.0)));
                        }

                    }
                }
            }
        }
    }


    molecule_free(&mol);

    return 0;
}

