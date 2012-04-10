#include <Stegosaurus.h>


__global__ void initDArrayKernel(double *m, int dim, double val) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim)
    m[idx] = val;
}

__global__ void compareMax(int dim, double *current_max, double *new_features) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    if (current_max[idx] < new_features[idx])
      current_max[idx] = new_features[idx];
  }
}

__global__ void compareMin(int dim, double *current_min, double *new_features) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    if (current_min[idx] > new_features[idx])
      current_min[idx] = new_features[idx];
  }
}

__global__ void rescaleKernel(int dim, double *vec_g, double *min_g, double *max_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim && max_g[idx] > 0.) {  // maybe need some invariant that deals with max=0 somehow, sets it to 1 for example
    vec_g[idx] = (vec_g[idx]-min_g[idx]) / max_g[idx];
  }
}

__global__ void varianceKernel(double divM, double* vec_g, double* mu_g, double* var_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  double delta = mu_g[idx] - vec_g[idx];
  
  var_g[idx] += delta * delta * divM;
}

__global__ void normalizeKernel(double *vec_g, double *mu_g, double *var_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  vec_g[idx] = (vec_g[idx] - mu_g[idx]) * var_g[idx];
}

inline void initDArray(double* m, int dim, int tpb, double val) {
  initDArrayKernel<<<BLOCKS(dim,tpb),tpb>>>(m, dim, val);
//   cudaThreadSynchronize();
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
  
  if (set->header->method == 0) {
    CUDA_CALL( cudaFree(set->max_g));
    CUDA_CALL( cudaFree(set->min_g));
    CUDA_CALL( cudaFree(set->mu_g));
    CUDA_CALL( cudaFree(set->var_g));
    CUDA_CALL( cudaFreeHost(&set->mu_vec));
  }
  if (set->counts != NULL)                CUDA_CALL( cudaFreeHost(set->counts));
  if (set->vec != NULL)                   CUDA_CALL( cudaFreeHost(set->vec));
  if (set->gauss->mu != NULL)             CUDA_CALL( cudaFree(set->gauss->mu));
  if (set->max_vec != NULL)               CUDA_CALL( cudaFreeHost(set->max_vec));
  if (set->qp_vec != NULL)                CUDA_CALL( cudaFreeHost(set->qp_vec));
  if (set->gauss->sigma != NULL)          free (set->gauss->sigma);
  if (set->gauss->sigma_inverse != NULL)  free (set->gauss->sigma_inverse);
  if (set->gauss->qr != NULL)             free (set->gauss->qr);
  if (set->gauss->qr_diag != NULL)        free (set->gauss->qr_diag);
  free(set->header);
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
  gp->doublesOnGPU = prop.totalGlobalMem * 3l / 4l / (long) sizeof(double);

  return gp;
}

stegoContext* init_stego() {
  stegoContext *steg = (stegoContext*) malloc(sizeof(stegoContext));
  steg->gpu_c = init_gpu();
  steg->features = NULL;
  steg->doublesInRAM = 5l * 1073741824l / (long) sizeof(double);

  return steg;
}


void close_stego(stegoContext *steg) {
  close_gpu(steg->gpu_c);
  closeFeatureSet(steg->features);
}

void close_gpu(gpuContext* gp) {
  cublasDestroy(gp->handle);
}

featureSet* openFeatureSet(const char *path, stegoContext *steg) {
//   printf("opening set \n");
  int i;
  uint64_t dim = 0ull;
  int read;
  uint64_t M;
  long posZero;
  int tpb = steg->gpu_c->threads_per_block;
//   double *vec;//, *qp, *max;// *mu_g, *vec_g, *max, *max_g, *qp_g, *qp;
  FILE *file = fopen(path, "r");
  featureHeader *header = (featureHeader*) malloc(sizeof(featureHeader));

  if (file == NULL) return NULL;
  featureSet *set = (featureSet*) malloc(sizeof(featureSet));
  readHeader(file, header);
  posZero = ftell(file);
  if (header->pair) {
    dim = (2*header->ranges[0][0]+1)*(2*header->ranges[0][1]+1) + 
          (2*header->ranges[1][0]+1)*(2*header->ranges[1][1]+1) + 
	  (2*header->ranges[2][0]+1)*(2*header->ranges[2][1]+1);
  } else {
    for (i = 0; i < 16; i++) dim += 2*header->ranges[0][i]+1;
    for (i = 0; i < 4; i++)  dim += 2*header->ranges[1][i]+1;
    for (i = 0; i < 15; i++) dim += 2*header->ranges[2][i]+1;
  }
  dim *= header->qp_range;
//   vec = (double*) malloc(dim*sizeof(double));
  
  set->header = header;
  set->dim = dim;
  set->files = (FILE**) malloc(MAX_FILES*sizeof(FILE*));//file;
  set->vsPerFile = (uint64_t*) malloc(MAX_FILES*sizeof(uint64_t));
  set->files[0] = file;
  set->num_files = 1;
  set->current_file = 0;
  set->dataOffset = posZero;
  set->gauss = (myGaussian*) malloc(sizeof(myGaussian));
  set->gauss->mu = NULL;
  set->gauss->sigma = NULL;
  set->gauss->sigma_inverse = NULL;
  set->gauss->qr = NULL;
  set->gauss->qr_diag = NULL;
  set->gauss->dim = dim;
  set->gpu_matrix_width = dim;                      // maybe wish to do something smarter here!
  CUDA_CALL( cudaHostAlloc(&set->counts, dim*sizeof(int), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->max_vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->min_vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->qp_vec, header->qp_range*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&set->ones_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&set->vec_g, dim*sizeof(double)));
  if (header->method == 0) {
    CUDA_CALL( cudaMalloc(&set->max_g, dim*sizeof(double)));
    CUDA_CALL( cudaMalloc(&set->min_g, dim*sizeof(double)));
    CUDA_CALL( cudaMalloc(&set->mu_g, dim*sizeof(double)));
    CUDA_CALL( cudaMalloc(&set->var_g, dim*sizeof(double)));
    CUDA_CALL( cudaHostAlloc(&set->mu_vec, dim*sizeof(double), cudaHostAllocDefault));
  }
  initDArray(set->ones_g, dim, tpb, 1.);
  
  M = 0ull;
  while ((read = fread(set->counts, sizeof(uint32_t), dim, file)) == dim) { // && M<10000
    M++;
  }
  if (read != 0) {
    printf("dim = &%i, read = %i \n", dim, read);
    printf("Wrong dimension?? \n");
  }

  set->vsPerFile[0] = M;
  set->M = M;
  set->divM = 1./(double) M;
  stegoRewind(set);
  
//   printf("esitmatig scaling \n");
  if (header->method == 0) {
    estimateScalingParameters(steg, set);
  }

//   printf("done. \n");
  return set;
}

int estimateScalingParameters(stegoContext *steg, featureSet *set) {
  uint64_t i;
  uint64_t M = set->M;
  uint64_t dim = set->dim;
  int tpb = steg->gpu_c->threads_per_block;
  
  initDArray(set->max_g, dim, tpb, 0.);
  initDArray(set->min_g, dim, tpb, INFINITY);
  initDArray(set->mu_g, dim, tpb, 0.);
  initDArray(set->var_g, dim, tpb, 0.);
  
//   printf("estimating scaling parameters over %ld vectors! \n", M);
  for (i = 0; i < M; i++) {
    readVectorL1D(steg, set, set->vec_g);
    compareMax<<<BLOCKS(dim,tpb),tpb>>>(dim, set->max_g, set->vec_g);
    compareMin<<<BLOCKS(dim,tpb),tpb>>>(dim, set->min_g, set->vec_g);
    cublasDaxpy(steg->gpu_c->handle, dim, &(set->divM), set->vec_g, 1, set->mu_g, 1);
  }
  stegoRewind(set);
  for (i = 0ull; i < M; i++) {
    readVectorL1D(steg, set, set->vec_g);
    varianceKernel<<<BLOCKS(dim,tpb),tpb>>>(set->divM, set->vec_g, set->mu_g, set->var_g);
  }
  stegoRewind(set);
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->max_g, 1, set->max_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->min_g, 1, set->min_vec, 1));
    CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->mu_g, 1, set->mu_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->var_g, 1, set->vec, 1));
  for (i = 0ull; i < dim; i++) {
    if (set->vec[i] > 0.)
      set->vec[i] = 1./sqrt(set->vec[i]);
  }
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), set->vec, 1, set->var_g, 1));
  for (i = 10000ull; i < 10100ull; i++) {
    printf("var[%i] = %g \n", i, set->vec[i]);
  }
}

int newFeatureFile(stegoContext *steg, featureSet* set, const char* path) {
  int i;
  int dim = set->dim;
  int read;
//   int vec[set->dim];
  uint64_t localM = 0ull;
  FILE *file;// = fopen(path, "r");
  featureHeader header;
  
//   printf("a new file, yay, paht is %s \n", path);
  if (set->num_files == MAX_FILES) return -1;
  file = fopen(path,"r");
  readHeader(file, &header);
//   printf("just read header, found method = %i\n", header.method);
  
  while ((read = fread(set->counts, sizeof(int), dim, file)) == dim) {// && M<10000
//     printf("1 vecotr :D, read %i elements \n", read);
    localM++;
  }
//   printf("it contains %ld vectors xD, read = %i \n", localM, read);
  set->M += localM;
  set->vsPerFile[set->num_files] = localM;
  set->divM = 1./set->M;
  set->files[set->num_files] = file;
  set->num_files++;
  stegoRewind(set);
  if (header.method == 0)
    estimateScalingParameters(steg, set);
  
  return 0;
}

// changes current file of set if necessary
int readCounts(featureSet *set) {
  int i;
  int read = 0;
  double *vec_g;
  double sum;
  
  read = fread(set->counts, sizeof(uint32_t), set->dim, set->files[set->current_file]);//readCountVector(vec, set->counts, set->dim, set->files[set->current_file]);
  if (read == 0) {
    set->current_file++;
    if (set->current_file == set->num_files)
      return -1;
    return readCounts(set);
  } else if (read != set->dim) {
    return -1;
  }
  
  for (i = 0; i < set->dim; i++) {
    set->vec[i] = (double) set->counts[i];
  }
  
  return read;
}

// // will do some uncompressing later, should only be called by other read___ methods!
// int readCountVector(double *data, int *cache, int dim, FILE *file) {
//   int i;
//   int read;
//   
//   read = fread(cache, sizeof(int), dim, file);
//   if (read == 0) return 0;
//   if (read != dim) return -1;
//   for (i = 0; i < read; i++) {
//     if (cache[i] < 0) printf("read something negative! \n");
//     data[i] = (double) cache[i];
//   }
//   return read;
// }

// reads directly into gpu memory
int readVectorL2(stegoContext *steg, featureSet *set, double *vec_g) {
  int read;
  double norm;
  
  read = readCounts(set);
  CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), set->vec, 1, vec_g, 1));
  CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, set->dim, vec_g, 1, &norm));
  norm = 1./norm;
  CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, set->dim, &norm, vec_g, 1));
//   CUBLAS_CALL( cublasGetVector(set->dim, sizeof(double), set->vec_g, 1, vec, 1));
  
  return read;
}

int readVectorL1D(stegoContext *steg, featureSet *set, double *vec_g) {
  int read;
  
  read = readCounts(set);
  scaleL1D(steg, set->dim, set->vec, vec_g, set->ones_g);
  
  if (read != set->dim) printf("wtf! \n");
  return read;
}

int readVectorRescaled(stegoContext *steg, featureSet *set, double *vec_g) {
  int read = readVectorL1D(steg, set, vec_g);
  int tpb = steg->gpu_c->threads_per_block;

//   CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), set->vec, 1, vec_g, 1));
//   scaleL1D(steg, set->dim, set->vec, vec_g, set->ones_g);
  rescaleKernel<<<BLOCKS(set->dim, tpb), tpb>>>(set->dim, vec_g, set->min_g, set->max_g);
  
  return read;
}

int readVectorNormalized(stegoContext *steg, featureSet *set, double *vec_g) {
  int read = readVectorL1D(steg, set, vec_g);
  int tpb = steg->gpu_c->threads_per_block;
  
  normalizeKernel<<<BLOCKS(set->dim, tpb), tpb>>>(vec_g, set->mu_g, set->var_g);
  
  return read;
}

// should write my own kernels and get rid of the ones_g
void scaleL1D(stegoContext *steg, int dim, double *vec, double *vec_g, double *ones_g) {
  double norm;
  
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), vec, 1, vec_g, 1));
  CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, dim, vec_g, 1, ones_g, 1, &norm));
  norm = (double) dim/norm;
  CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, dim, &norm, vec_g, 1));
}
