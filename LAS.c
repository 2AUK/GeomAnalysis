#include <math.h>

#include "LAS.h"
#include "molecule.h"

double unit_vector(Molecule mol,int cart, int a, int b){
    return -1 * ((mol.geom[a + cart] - mol.geom[b + cart]) / molecule_bond(a, b)); 
}
