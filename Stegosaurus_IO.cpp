#include <Stegosaurus.h>

featureSet* openFeatureSet(char *path){
  int dim = 0;
  int read;
  int M;
  double *vec;
  FILE *file = fopen(path, "r");
  featureSet *set = (featureSet*) malloc(sizeof(featureSet));
  
  // scan throgh file once to get number of features
  fread(&dim, sizeof(int), 1, file);
  vec = (double*) malloc(dim*sizeof(double));
  M = 0;
  while ((read = fread(vec, sizeof(double), dim, file))>0)
    M++;
  if (read == -1)
    return NULL;
  
  set->dim = dim;
  set->M = M;
  set->file = file;
  set->current_feature = vec;
  
  stegoRewind(set);
//   printf("opened feature set: dim=%i \n", dim);
  
  return set;
}

int closeFeatureSet(featureSet *set) {
  int i;
  i = fclose(set->file);
  free(set->current_feature);
  free(set);
  return i;
}

int readVector(featureSet *set, double *vec) {
  int read;
  
  read = fread(vec, sizeof(double), set->dim, set->file);
//   for (int i = 0; i < 2560; i++) {
//     printf("%f ", vec[i]);
//   }
//   printf("\n");
  if (read > 0 && read != set->dim)
    return -1;
  return read;
}

void stegoRewind(featureSet *set) {
  int dim;
  
  rewind(set->file);
  fread(&dim, sizeof(int), 1, set->file);
}

void storeGaussian(char* path, myGaussian* gauss) {
  
}
