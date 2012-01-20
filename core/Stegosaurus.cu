#include <Stegosaurus.h>


featureSet* openFeatureSet(char *path) {
  int dim = 0;
  int read;
  int M;
  double *vec, *qp, *max;// *mu_g, *vec_g, *max, *max_g, *qp_g, *qp;
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
  CUDA_CALL( cudaHostAlloc(&max, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&qp, QP_RANGE*sizeof(double), cudaHostAllocDefault));
  vec = (double*) malloc(dim*sizeof(double));
  
  M = 0;
  while ((read = fread(vec, sizeof(double), dim, file))>0) // && M<10000
    M++;
  if (read == -1)
    return NULL;
  
  set->dim = dim;
  set->M = M;
  set->divM = 1./((double) M);
  set->kl_div = -1.;                                // some invalid value
  set->file = file;
//   set->current_feature = vec;
  set->max_vec = max;
  set->qp_vec = qp;
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
  set->gpu_matrix_width = dim;                      // maybe wish to do something smarter here!
  
  stegoRewind(set);
//   printf("opened feature set: dim=%i, M=%i \n", dim, M);
  free(vec);
  
  return set;
}

int closeFeatureSet(featureSet *set) {
  int i;
  i = fclose(set->file);

//   CUDA_CALL( cudaFreeHost(set->current_feature));
  CUDA_CALL( cudaFreeHost(set->max_vec));
  CUDA_CALL( cudaFreeHost(set->qp_vec));
//   CUDA_CALL( cudaFree(set->mu_g));
//   CUDA_CALL( cudaFree(set->qp_g));
//   CUDA_CALL( cudaFree(set->vec_g));
//   CUDA_CALL( cudaFree(set->max_g));
  if (set->gauss->mu != NULL)             free (set->gauss->mu);
  if (set->gauss->sigma != NULL)          free (set->gauss->sigma);
  if (set->gauss->sigma_inverse != NULL)  free (set->gauss->sigma_inverse);
  if (set->gauss->qr != NULL)             free (set->gauss->qr);
  if (set->gauss->qr_diag != NULL)        free (set->gauss->qr_diag);
  free(set->gauss);
  free(set->name);
  free(set);
  return i;
}

gpuContext* init_gpu() {
  cudaDeviceProp prop;
  gpuContext *gp = (gpuContext*) malloc(sizeof(gpuContext));

  CUDA_CALL( cudaGetDeviceProperties(&prop, 0));
  CUBLAS_CALL( cublasCreate(&(gp->handle)));
  gp->threads_per_block = prop.maxThreadsPerBlock;
//   printf("%i \n", gp->threads_per_block);

  return gp;
}

stegoContext* init_stego() {
  stegoContext *steg = (stegoContext*) malloc(sizeof(stegoContext));
  steg->gpu_c = init_gpu();
  steg->features = NULL;

  return steg;
}


void close_stego(stegoContext *steg) {
  close_gpu(steg->gpu_c);
  closeFeatureSet(steg->features);
}

void close_gpu(gpuContext* gp) {
  cublasDestroy(gp->handle);
}
