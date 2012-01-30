#include "Stegosaurus.h"


featureSet* openFeatureSet(char *path) {
  int dim = 0;
  int read;
  int M;
  double *vec;//, *qp, *max;// *mu_g, *vec_g, *max, *max_g, *qp_g, *qp;
  FILE *file = fopen(path, "r");
//   cudaPointerAttributes attr;
  
//   printf("trying to open: %s \n", path);
  if (file == NULL) return NULL;
  fread(&dim, sizeof(int), 1, file);
  
//   printf("Opening some not NULL featureset!, dim=%i \n", dim);
  featureSet *set = (featureSet*) malloc(sizeof(featureSet));
//   CUDA_CALL( cudaMalloc(&mu_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&qp_g, 20*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&max_g, dim*sizeof(double)));
//   CUDA_CALL( cudaHostAlloc(&vec, dim*sizeof(double), cudaHostAllocDefault));
//   CUDA_CALL( cudaHostAlloc(&max, dim*sizeof(double), cudaHostAllocDefault));
//   CUDA_CALL( cudaHostAlloc(&qp, QP_RANGE*sizeof(double), cudaHostAllocDefault));
  vec = (double*) malloc(dim*sizeof(double));
  
  M = 0;
  while ((read = fread(vec, sizeof(double), dim, file))>0) // && M<10000
    M++;
//   if (read == -1)
//     return NULL;
  
  set->dim = dim;
  set->M = M;
  set->divM = 1./((double) M);
  set->kl_div = -1.;                                // some invalid value
  set->files = (FILE**) malloc(MAX_FILES*sizeof(FILE*));//file;
  set->files[0] = file;
  set->num_files = 1;
  set->current_file = 0;
//   set->current_feature = vec;
//   set->max_vec = max;
//   set->qp_vec = qp;
//   set->mu_g = mu_g;
//   set->qp_g = qp_g;
//   set->vec_g = vec_g;
//   set->max_g = max_g;
  set->gauss = (myGaussian*) malloc(sizeof(myGaussian));
  set->gauss->mu = NULL;
  set->gauss->sigma = NULL;
  set->gauss->sigma_inverse = NULL;
  set->gauss->qr = NULL;
  set->gauss->qr_diag = NULL;
  set->gauss->dim = dim;
  set->gpu_matrix_width = dim/2+1;                      // maybe wish to do something smarter here!
  
  stegoRewind(set);
//   printf("opened feature set: dim=%i, M=%i \n", dim, M);
  free(vec);
  
  return set;
}

int addFeatureFile(featureSet *set, char *path) {
  int dim;
  int read;
  double vec[set->dim];
  FILE *file;
  
  printf("adding featuer file \n");
  
  if (set->num_files == MAX_FILES) return -1;
  file = fopen(path,"r");
  
  fread(&dim, sizeof(int), 1, file);
  if (dim != set->dim) {
    printf("Dimension mismatch! \n");
    fclose(file);
    return -2;
  }
  while ((read = fread(vec, sizeof(double), dim, file))>0) // && M<10000
    set->M++;
  set->divM = 1./set->M;
  rewind(file);
  
  set->files[set->num_files] = file;
  set->num_files++;
  return 0;
}

int readVector(featureSet *set, double *vec) {
  int read = 0;
  int dim;
  
  while (read == 0 && set->current_file < set->num_files) {
    read = fread(vec, sizeof(double), set->dim, set->files[set->current_file]);
    if (read > 0 && read != set->dim)
      return -1;
    else if (read == 0) {
      set->current_file++;
      fread(&dim, sizeof(int), 1, set->files[set->current_file]);
    }
  }
  return read;
}

void stegoRewind(featureSet *set) {
  int i;
  int dim;
  
  for (i = 0; i < set->num_files; i++) {
    rewind(set->files[i]);
  }
  set->current_file = 0;
  fread(&dim, sizeof(int), 1, set->files[set->current_file]);
}

// void storeGaussian(char* path, myGaussian* gauss) {
//   
// }
