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

__global__ void rescaleVec(int dim, double *vec_g, double *min_g, double *max_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim && max_g[idx] > 0.) {  // maybe need some invariant that deals with max=0 somehow, sets it to 1 for example
    vec_g[idx] = (vec_g[idx]-min_g[idx]) / max_g[idx];
  }
}

inline void initDArray(double* m, int dim, int tpb, double val) {
  initDArrayKernel<<<BLOCKS(dim,tpb),tpb>>>(m, dim, val);
  cudaThreadSynchronize();
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
  
  if (set->counts != NULL)                CUDA_CALL( cudaFreeHost(set->counts));
  if (set->vec != NULL)                   CUDA_CALL( cudaFreeHost(set->vec));
  if (set->gauss->mu != NULL)             CUDA_CALL( cudaFree(set->gauss->mu));
  if (set->max_vec != NULL)               CUDA_CALL( cudaFreeHost(set->max_vec));
  if (set->qp_vec != NULL)                CUDA_CALL( cudaFreeHost(set->qp_vec));
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
  int i;
  int dim = 0;
  int read;
  int M;
  long posZero;
  int tpb = steg->gpu_c->threads_per_block;
//   double *vec;//, *qp, *max;// *mu_g, *vec_g, *max, *max_g, *qp_g, *qp;
  FILE *file = fopen(path, "r");
  featureHeader header;

  if (file == NULL) return NULL;
  featureSet *set = (featureSet*) malloc(sizeof(featureSet));
  readHeader(file, &header);
  posZero = ftell(file);
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
//   vec = (double*) malloc(dim*sizeof(double));
  
  set->dim = dim;
  set->files = (FILE**) malloc(MAX_FILES*sizeof(FILE*));//file;
  set->vsPerFile = (long*) malloc(MAX_FILES*sizeof(long));
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
  set->gpu_matrix_width = dim/2+1;                      // maybe wish to do something smarter here!
  CUDA_CALL( cudaHostAlloc(&set->counts, dim*sizeof(int), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->max_vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->min_vec, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&set->qp_vec, header.qp_range*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&set->ones_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&set->max_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&set->min_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&set->vec_g, dim*sizeof(double)));
  initDArray(set->max_g, dim, tpb, 0.);
  initDArray(set->min_g, dim, tpb, INFINITY);
  initDArray(set->ones_g, dim, tpb, 1.);
  
  M = 0;
  while ((read = readCountVector(set->vec, set->counts, dim, file))>0) { // && M<10000
    // maybe want to find min or max here?
    scaleL1D(steg, dim, set->vec, set->vec_g, set->ones_g);
    compareMax<<<BLOCKS(dim,tpb),tpb>>>(dim, set->max_g, set->vec_g);
    compareMin<<<BLOCKS(dim,tpb),tpb>>>(dim, set->min_g, set->vec_g);
    M++;
  }
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->max_g, 1, set->max_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->min_g, 1, set->min_vec, 1));
  if (read != 0) printf("Wrong dimension?? \n");
  rewind(file);                           // need a method that
  readHeader(file, &header);              // does this using fseek!
  
  set->M = M;
  set->divM = 1./(double) M;
  
//   free(vec);
  
  return set;
}

// changes current file of set if necessary
int readCounts(stegoContext *steg, featureSet *set, double *vec) {
  int i;
  int read = 0;
  double *vec_g;
  double sum;
  
  while (read == 0 && set->current_file < set->num_files) {
    read = readCountVector(vec, set->counts, set->dim, set->files[set->current_file]);
//     if (read != set->dim) printf("Read the following number of bytes: %i \n", read);
    
    if (read != set->dim)
      return -1;
    else if (read == 0) {
      set->current_file++;
    }
  }
  
  return read;
}

// will do some uncompressing later, should only be called by other read___ methods!
int readCountVector(double *data, int *cache, int dim, FILE *file) {
  int i;
  int read;
  
  read = fread(cache, sizeof(int), dim, file);
  if (read == 0) return 0;
  if (read != dim) return -1;
  for (i = 0; i < read; i++) {
    data[i] = (double) cache[i];
  }
  return read;
}

// reads directly into gpu memory
int readVectorL2(stegoContext *steg, featureSet *set, double *vec_g) {
  int read;
  double norm;
  
  read = readCounts(steg, set, set->vec);
  CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), set->vec, 1, vec_g, 1));
  CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, set->dim, vec_g, 1, &norm));
  norm = 1./norm;
  CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, set->dim, &norm, vec_g, 1));
//   CUBLAS_CALL( cublasGetVector(set->dim, sizeof(double), set->vec_g, 1, vec, 1));
  
  return read;
}

int readVectorL1D(stegoContext *steg, featureSet *set, double *vec_g) {
  int read;
//   double norm;
  
  read = readCounts(steg, set, set->vec);
  scaleL1D(steg, set->dim, set->vec, vec_g, set->ones_g);
//   CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), set->vec, 1, vec_g, 1));
//   CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, set->dim, vec_g, 1, set->ones_g, 1, &norm));
//   norm = (double) set->dim/norm;
//   CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, set->dim, &norm, vec_g, 1));
  
  return read;
}

int readVectorRescaled(stegoContext *steg, featureSet *set, double *vec_g) {
  int read = readCounts(steg, set, set->vec);
  int tpb = steg->gpu_c->threads_per_block;

  CUBLAS_CALL( cublasSetVector(set->dim, sizeof(double), set->vec, 1, vec_g, 1));
  scaleL1D(steg, set->dim, set->vec, vec_g, set->ones_g);
  rescaleVec<<<BLOCKS(set->dim, tpb), tpb>>>(set->dim, vec_g, set->min_g, set->max_g);
  
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
