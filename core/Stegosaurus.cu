#include <Stegosaurus.h>

__global__ void initDArray2(double *m, int dim, double val) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim)
    m[idx] = val;
}

void printPointerInfo(void *ptr) {
  cudaPointerAttributes attr;
  
  cudaPointerGetAttributes(&attr, ptr);
  if (attr.memoryType == cudaMemoryTypeHost) printf("Pointer type: Host \n");
  else if (attr.memoryType == cudaMemoryTypeDevice) printf("Pointer type: Device \n");
  else printf("Pointer Attr: ?? \n");
}

int closeFeatureSet(featureSet *set) {
  int i;
  
  for (i = 0; i < set->num_files; i++) {
    fclose(set->files[i]);
  }

//   CUDA_CALL( cudaFreeHost(set->current_feature));
//   CUDA_CALL( cudaFreeHost(set->max_vec));
//   CUDA_CALL( cudaFreeHost(set->qp_vec));
//   CUDA_CALL( cudaFree(set->mu_g));
//   CUDA_CALL( cudaFree(set->qp_g));
//   CUDA_CALL( cudaFree(set->vec_g));
//   CUDA_CALL( cudaFree(set->max_g));
//   if (set->gauss->mu != NULL)             free (set->gauss->mu);
  if (set->gauss->mu != NULL)             CUDA_CALL( cudaFree(set->gauss->mu));
  if (set->max_vec != NULL)               CUDA_CALL( cudaFree(set->max_vec));
  if (set->qp_vec != NULL)                CUDA_CALL( cudaFree(set->qp_vec));
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
  gp->num_streams = 4;
//   gp->qr_cache = 10;
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
//   printf("Opening featureset with dim=%i \n", dim);
  
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
  while ((read = fread(vec, sizeof(int), dim, file))>0) // && M<10000
    M++;
//   printf("M = %i \n", M);
  if (read != 0) printf("Wrong dimension?? \n");
  rewind(file);
  readHeader(file, &header);
//   if (read == -1)
//     return NULL;
  
  set->dim = dim;
  set->M = M;
  set->divM = 1./((double) M);
//   set->kl_div = -1.;                                // some invalid value
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
  CUDA_CALL( cudaMalloc(&set->vec_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&set->ones_g, dim*sizeof(double)));
  initDArray2<<<BLOCKS(dim, 1024), 1024>>>(set->ones_g, dim, 1.);
  
//   stegoRewind(set);
//   printf("opened feature set: dim=%i, M=%i \n", dim, M);
  free(vec);
  
  return set;
}

int readVector(stegoContext *steg, featureSet *set, double *vec) {
  int i;
  int read = 0;
//   featureSet *set = steg->features;
  int vec_i[set->dim];
  double *vec_g;
  double sum;
  
  while (read == 0 && set->current_file < set->num_files) {
    read = fread(vec_i, sizeof(int), set->dim, set->files[set->current_file]);
    for (i = 0; i < set->dim; i++) {
      vec[i] = (double) vec_i[i];
    }
    CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), vec, 1, set->vec_g, 1));
    CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, set->dim, set->vec_g, 1, set->ones_g, 1, &sum));
//     printf("sum = %g \n", sum);
    sum = ((double) set->dim)/sum;
    CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, set->dim, &sum, set->vec_g, 1));
    CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, set->dim, set->vec_g, 1, set->ones_g, 1, &sum));
//     printf("sum = %g \n", sum);
    CUBLAS_CALL( cublasGetVector(set->dim, sizeof(double), set->vec_g, 1, vec, 1));
    
    if (read > 0 && read != set->dim)
      return -1;
    else if (read == 0) {
      set->current_file++;
//       fread(&dim, sizeof(int), 1, set->files[set->current_file]);
    }
  }
  
//   CUDA_CALL();
  
  return read;
}

