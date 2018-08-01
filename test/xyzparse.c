#include <stdio.h>
#include <stdlib.h>

#define BUF_SIZE 10000

void readXYZ(FILE *stream, int *z_array, double *coord_array, int *natom);
void printZAr(int *z_array, int size);

int main(int argc, char *argv[]){
    int *z_array;
    double *coord_array;
    int natom;
    FILE* input_molecule;

    input_molecule = fopen("acetaldehyde.txt", "r");
    readXYZ(input_molecule, z_array, coord_array, &natom);
    printf("%d\n", natom);
    printZAr(z_array, natom);
    return 0;
}

void readXYZ(FILE *stream, int *z_array, double *coord_array, int *natom){
    char c;
    char cnatom[3];
    int count, i, j;
    count = i = j = 0;
    double z_val,x,y,z;
    char *buf;

    while ((c = fgetc(stream)) != '\n')
        cnatom[count++] = c;
    *natom = atoi(cnatom);

   coord_array = (double*) malloc(sizeof(double) * *natom * *natom);
   z_array = (int*) malloc(sizeof(int) * *natom);
   buf = (char*) malloc(sizeof(char) * 4 * *natom);

   while (fgets(buf, BUF_SIZE, stream) != NULL){
       sscanf(buf, "%lf %lf %lf %lf", &z_val, &x, &y, &z);
       z_array[i++] = z_val;
       coord_array[j++] = x;
       coord_array[j++] = y;
       coord_array[j++] = z;
   }
   printf("%d\n", *natom);
}

void printZAr(int *z_array, int size){
    for (int i = 0; i < size; i++)
       printf("%d\t", z_array[i]);
   putchar('\n'); 
}
