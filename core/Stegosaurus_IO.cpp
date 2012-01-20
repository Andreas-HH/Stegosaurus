#include "Stegosaurus.h"


int readVector(featureSet *set, double *vec) {
  int read;
  
  read = fread(vec, sizeof(double), set->dim, set->file);
  if (read > 0 && read != set->dim)
    return -1;
  return read;
}

void stegoRewind(featureSet *set) {
  int dim;
  
  rewind(set->file);
  fread(&dim, sizeof(int), 1, set->file);
}

// void storeGaussian(char* path, myGaussian* gauss) {
//   
// }
