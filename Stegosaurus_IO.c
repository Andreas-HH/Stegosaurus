#include<Stegosaurus_IO.h>

FeatureSet* openFeatureSet(char *path){
  int dim;
  double *vec;
  FILE *file = fopen(path, "r");
  FeatureSet *set = malloc(sizeof(FeatureSet));
  
  dim = fread(&dim, sizeof(int), 1, file);
  vec = malloc(dim*sizeof(double));
  set->dim = dim;
  set->file = file;
  set->current_feature = vec;
  
  return set;
}

int closeFeatureSet(FeatureSet *set) {
  int i;
  i = fclose(set->file);
  free(set->current_feature);
  free(set);
  return i;
}

int readVector(FeatureSet *set) {
  int read;
  
  read = fread(set->current_feature, sizeof(double), set->dim, set->file);
  if (read > 0 && read != set->dim)
    return -1;
  return read;
}

void storeGaussian(char* path, MyGaussian* gauss) {
  
}
