#include<Stegosaurus_Estimate.h>

int estimateGaussian(char *features) {
  int i;
  int read;
  int M; // # of vectors
  
  init_gpu();
  FeatureSet *set = openFeatureSet(features);
  
  M = 0;
  while ((read = readVector(set))>0) {
    M++;
  }
  if (read == -1)
    return -1;
  rewind(set->file);
  printf("%i, %i \n", read, M);
  
  closeFeatureSet(set);
}