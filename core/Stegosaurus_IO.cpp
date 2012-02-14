#include "Stegosaurus.h"


int readHeader(FILE *file, featureHeader *header) {
  int i;
  
  printf("pair: %i \n", header->pair);
  fread(&header->pair, sizeof(char), 1, file);
  printf("pair: %i \n", header->pair);
  fread(&header->slice_type, sizeof(char), 1, file);
  printf("slice_type: %i \n", header->slice_type);
  fread(&header->method, sizeof(char), 1, file);
  printf("method: %i \n", header->method);
  if (header->method != 0) {
    fread(&header->using_rate, sizeof(char), 1, file);
    fread(&header->rate, sizeof(double), 1, file);
    printf("rate: %f \n", header->rate);
    fread(&header->accept, sizeof(char), 1, file);
  }
  fread(&header->qp_offset, sizeof(char), 1, file);
  printf("offset: %i \n", header->qp_offset);
  fread(&header->qp_range, sizeof(char), 1, file);
  printf("range: %i \n", header->qp_range);
  
  for (i = 0; i < 16; i++) fread(&header->ranges[0][0]+i, sizeof(unsigned char), 1, file);
  for (i = 0; i < 4; i++)  fread(&header->ranges[1][0]+i, sizeof(unsigned char), 1, file);
  for (i = 0; i < 15; i++) fread(&header->ranges[2][0]+i, sizeof(unsigned char), 1, file);
  
  return 0;
}

featureSet* openFeatureSet(const char *path) {
  int i;
  int dim = 0;
  int read;
  int M;
  double *vec;//, *qp, *max;// *mu_g, *vec_g, *max, *max_g, *qp_g, *qp;
  FILE *file = fopen(path, "r");
  featureHeader header;
//   cudaPointerAttributes attr;
  
//   printf("trying to open: %s \n", path);
  if (file == NULL) return NULL;
  featureSet *set = (featureSet*) malloc(sizeof(featureSet));
  readHeader(file, &header);
  if (header.pair) {
    dim = (2*header.ranges[0][0]+1)*(2*header.ranges[0][1]+1) + 
          (2*header.ranges[1][0]+1)*(2*header.ranges[1][1]+1) + 
	  (2*header.ranges[2][0]+1)*(2*header.ranges[2][1]+1);
  } else {
    for (i = 0; i < 16; i++) dim += 2*header.ranges[0][i]+1;
    for (i = 0; i < 4; i++)  dim += 2*header.ranges[0][i]+1;
    for (i = 0; i < 15; i++) dim += 2*header.ranges[0][i]+1;
  }
  dim *= header.qp_range;
  printf("Opening featureset with dim=%i \n", dim);
  
//   printf("Opening some not NULL featureset!, dim=%i \n", dim);
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
  printf("M = %i \n", M);
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

int newFeatureFile(featureSet* set, const char* path) {
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
//       fread(&dim, sizeof(int), 1, set->files[set->current_file]);
    }
  }
  return read;
}

void stegoRewind(featureSet *set) {
  int i;
  featureHeader header;
  
  for (i = 0; i < set->num_files; i++) {
    rewind(set->files[i]);
    readHeader(set->files[i], &header);
  }
  set->current_file = 0;
//   fread(&dim, sizeof(int), 1, set->files[set->current_file]);
}

// void storeGaussian(char* path, myGaussian* gauss) {
//   
// }
