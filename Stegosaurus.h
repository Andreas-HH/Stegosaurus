#ifndef STEGOSAURUS
#define STEGOSAURUS

#include <stdio.h>
#include <stdlib.h>

typedef struct FeatureSet {
  int dim;
  FILE *file;
  double *current_feature;
} FeatureSet;

typedef struct MyGaussian {
//   int hasQR;
//   int hasInverse;
  int dim;
  double *mu;
  double *sigma;
  double *sigma_inverse;
  double *qr_diag;
  double *qr;
//   int dim;
} MyGaussian;

typedef struct stegoContext {
  
} stegoContext;

#endif /* STEGOSAURUS */