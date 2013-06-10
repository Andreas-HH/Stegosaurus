#include "Stegosaurus.h"
#include "StegoClassifier.h"

__global__ void initDArrayKernel(double *m, int dim, double val) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim)
    m[idx] = val;
}

__global__ void finishMax(int dim, double *min, double *max) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    max[idx] = max[idx] - min[idx];
    if (max[idx] < 0.0000001)
      max[idx] = 1.;
  }
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

// same as normalizeKernel? even safer
__global__ void rescaleKernel(int dim, double *vec_g, double *min_g, double *max_g) { 
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim && max_g[idx] > 0.) {  // maybe need some invariant that deals with max=0 somehow, sets it to 1 for example
    vec_g[idx] = (vec_g[idx]-min_g[idx]) / max_g[idx];
  }
}

__global__ void varianceKernel(double divM, double* vec_g, double* mu_g, double* var_g, int dim) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  double delta = mu_g[idx] - vec_g[idx];
  
  if (idx > dim)
    var_g[idx] += delta * delta * divM;
}

__global__ void normalizeKernel(double *vec_g, double *mu_g, double *var_g, int dim) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim)
    vec_g[idx] = (vec_g[idx] - mu_g[idx]) * var_g[idx];
}

void initDArray(double* m, int dim, int tpb, double val) {
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
    if (set->files[i] != 0)
      fclose(set->files[i]);
  }
  free(set->files);
  free(set->paths);
  
  if (set->header->method == 0) {
    CUDA_CALL( cudaFree(set->max_g));
    CUDA_CALL( cudaFree(set->min_g));
    CUDA_CALL( cudaFree(set->mu_g));
    CUDA_CALL( cudaFree(set->var_g));
    CUDA_CALL( cudaFreeHost(&set->mu_vec));
    CUDA_CALL( cudaFreeHost(&set->mask_vec));
    CUDA_CALL( cudaFreeHost(&set->mask_counts));
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
  int dim = 0, dim_file = 0;
  int hist_dim = 0;
  int pair_dim = 0;
  int uvsv_dim = 0;
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
  for (i = 0; i < 16; i++) hist_dim += 2*header->ranges[0][i]+1;
  for (i = 0; i < 4; i++)  hist_dim += 2*header->ranges[1][i]+1;
  for (i = 0; i < 15; i++) hist_dim += 2*header->ranges[2][i]+1;
  pair_dim = (2*header->ranges[0][0]+1)*(2*header->ranges[0][1]+1) + 
             (2*header->ranges[1][0]+1)*(2*header->ranges[1][1]+1) + 
	     (2*header->ranges[2][0]+1)*(2*header->ranges[2][1]+1);
  uvsv_dim += (2*header->ranges[1][0]+1)*(2*header->ranges[1][0]+1);
  for (i = 1; i < 16; i++) {
    if (2*header->ranges[2][i-1] > 0) uvsv_dim += (2*header->ranges[2][i-1]+1)*(2*header->ranges[2][i-1]+1);
  }
//   if (header->pair) {
//     dim = (2*header->ranges[0][0]+1)*(2*header->ranges[0][1]+1) + 
//           (2*header->ranges[1][0]+1)*(2*header->ranges[1][1]+1) + 
// 	  (2*header->ranges[2][0]+1)*(2*header->ranges[2][1]+1);
//   } else {
//     for (i = 0; i < 16; i++) dim += 2*header->ranges[0][i]+1;
//     for (i = 0; i < 4; i++)  dim += 2*header->ranges[1][i]+1;
//     for (i = 0; i < 15; i++) dim += 2*header->ranges[2][i]+1;
//   }
  dim = (int)header->qp_range * (hist_dim + pair_dim + uvsv_dim);
  dim_file = (int)header->qp_range * (hist_dim + pair_dim + uvsv_dim);
//   printf("have dim: %i, %i, %i \n", hist_dim, pair_dim, uvsv_dim);
//   vec = (double*) malloc(dim*sizeof(double));
  
  set->header = header;
  set->hist_dim = hist_dim;
  set->pair_dim = pair_dim;
  set->uvsv_dim = uvsv_dim;
  set->dim = dim;
  set->dim_file = dim_file;
  set->files = (FILE**) malloc(MAX_FILES*sizeof(FILE*));//file;
  set->paths = (char**) malloc(MAX_FILES*sizeof(char*));
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
  CUDA_CALL( cudaHostAlloc(&set->counts, dim_file*sizeof(store_elem), cudaHostAllocDefault));
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
    CUDA_CALL( cudaHostAlloc(&set->mask_counts, dim*sizeof(uint64_t), cudaHostAllocDefault));
    CUDA_CALL( cudaHostAlloc(&set->mask_vec, dim*sizeof(int), cudaHostAllocDefault));
  }
  initDArray(set->ones_g, dim, tpb, 1.);
  
  M = 0ull;
  while ((read = fread(set->counts, sizeof(store_elem), dim_file, file)) == dim_file) { // && M<10000
    M++;
  }
  if (read != 0) {
    printf("dim = %d, read = %i \n", dim, read);
    printf("Wrong dimension?? \n");
  }

  set->paths[0] = (char*) malloc((strlen(path)+1)*sizeof(char));
  memcpy(set->paths[0], path, (strlen(path)+1)*sizeof(char));
  set->vsPerFile[0] = M;
  set->M = M;
  set->divM = 1./(double) M;
  stegoRewind(set);
  
//   printf("esitmatig scaling \n");
  if (header->method == 0) {
    estimateScalingParameters(steg, set);
  }
  fclose(file);
  set->files[0] = 0;

//   printf("done. \n");
  return set;
}

int estimateScalingParameters(stegoContext *steg, featureSet *set) {
  uint64_t i, j;
  uint64_t M = set->M;
  uint64_t dim = set->dim;
  uint64_t max_elem = 0ul;
  int tpb = steg->gpu_c->threads_per_block;
  
  initDArray(set->max_g, dim, tpb, 0.);
  initDArray(set->min_g, dim, tpb, INFINITY);
  initDArray(set->mu_g, dim, tpb, 0.);
  initDArray(set->var_g, dim, tpb, 0.);
  for (j = 0ull; j < set->dim; j++) {
    set->mask_counts[j] = 0ul;
    set->mask_vec[j] = 0;
  }
  
  for (i = 0ull; i < M; i++) {
    readVectorL1D(steg, set, set->vec_g);
    for (j = 0ull; j < set->dim; j++) {
      if (set->counts[j] > 0ul)
	set->mask_counts[j]++;
      if (set->counts[j] > max_elem) 
	max_elem = set->counts[j];
    }
    compareMax<<<BLOCKS(dim,tpb),tpb>>>(dim, set->max_g, set->vec_g);
    compareMin<<<BLOCKS(dim,tpb),tpb>>>(dim, set->min_g, set->vec_g);
    cublasDaxpy(steg->gpu_c->handle, dim, &(set->divM), set->vec_g, 1, set->mu_g, 1);
  }
  finishMax<<<BLOCKS(dim,tpb),tpb>>>(dim, set->min_g, set->max_g);
  for (j = 0ull; j < set->dim; j++) {
    if (set->mask_vec[j] > set->M/100ull) set->mask_vec[j] = 1;
  }
  printf("max_elem: %u \n", max_elem);
  stegoRewind(set);
  for (i = 0ull; i < M; i++) {
    readVectorL1D(steg, set, set->vec_g);
    varianceKernel<<<BLOCKS(dim,tpb),tpb>>>(set->divM, set->vec_g, set->mu_g, set->var_g, dim);
  }
  stegoRewind(set);
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->max_g, 1, set->max_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->min_g, 1, set->min_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->mu_g, 1, set->mu_vec, 1));
  // maybe this can be done in a kernel!!
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), set->var_g, 1, set->vec, 1));
  for (i = 0ull; i < dim; i++) {
    if (set->vec[i] > 0.)
      set->vec[i] = 1./sqrt(set->vec[i]);
    else {
      set->vec[i] = 1.;
    }
  }
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), set->vec, 1, set->var_g, 1));

  return 0;
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
  
  while ((read = fread(set->counts, sizeof(store_elem), set->dim_file, file)) == set->dim_file) {// && M<10000
//     printf("1 vecotr :D, read %i elements \n", read);
    localM++;
  }
  fclose(file);
//   printf("it contains %ld vectors xD, read = %i \n", localM, read);
  set->M += localM;
  set->vsPerFile[set->num_files] = localM;
  set->divM = 1./set->M;
  set->files[set->num_files] = 0;
  set->paths[set->num_files] = (char*) malloc((strlen(path)+1)*sizeof(char));
  memcpy(set->paths[set->num_files], path, (strlen(path)+1)*sizeof(char));
  set->num_files++;
//   stegoRewind(set);
  if (header.method == 0) {
    startAction(set);
    estimateScalingParameters(steg, set);
    endAction(set);
  }
  
//   for (i = 0; i < set->num_files; i++) {
//     printf("\"%s\" (%s), %i \n", set->paths[i], path, strlen(path));
//   }
  
  return 0;
}

// changes current file of set if necessary
int readCounts(featureSet *set) {
  int i;
  int read = 0;
  
  read = fread(set->counts, sizeof(store_elem), set->dim_file, set->files[set->current_file]);//readCountVector(vec, set->counts, set->dim, set->files[set->current_file]);
  if (read == 0) {
    fseek(set->files[set->current_file], set->dataOffset, SEEK_SET);
    set->current_file++;
    if (set->current_file == set->num_files)
      return -1;
    fseek(set->files[set->current_file], set->dataOffset, SEEK_SET);
    return readCounts(set);
  } else if (read != set->dim_file) {
    return -1;
  }
  
  for (i = 0; i < set->dim; i++) {
    set->vec[i] = (double) set->counts[i];
    if (set->vec[i] < 0.) printf(";_; \n");
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
  
  if (read != set->dim_file) printf("read something wrong! %i, %i \n", set->dim_file, read);
  return read;
}

int readVectorRescaled(stegoContext *steg, featureSet *set, double *vec_g) {
  int read = readVectorL1D(steg, set, vec_g);
  int tpb = steg->gpu_c->threads_per_block;

  rescaleKernel<<<BLOCKS(set->dim, tpb), tpb>>>(set->dim, vec_g, set->min_g, set->max_g);
  
  return read;
}

int readVectorNormalized(stegoContext *steg, featureSet *set, double *vec_g) {
  int read = readVectorL1D(steg, set, vec_g);
  int tpb = steg->gpu_c->threads_per_block;
  
  normalizeKernel<<<BLOCKS(set->dim, tpb), tpb>>>(vec_g, set->mu_g, set->var_g, set->dim);
  
  return read;
}

void scaleL1D(stegoContext *steg, int dim, double *vec, double *vec_g, double *ones_g) {
  double norm;
  
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), vec, 1, vec_g, 1));
  CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, dim, vec_g, 1, ones_g, 1, &norm));
  norm = (double) dim/norm;
  CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, dim, &norm, vec_g, 1));
}





// #include "Stegosaurus.h"


void pathConcat(const char* a, const char* b, char* result) {
  int i;
  int count;
  
  for (i = 0; a[i] != '\0'; i++)
    result[i] = a[i];
  result[i] = '/';
  count = i+1;
  for (i = 0; b[i] != '\0'; i++)
    result[count+i] = b[i];
  result[count+i] = '\0';
}


FeatureCollection::FeatureCollection(featureHeader *h) {
  memcpy(&header, h, sizeof(featureHeader));
}

FeatureCollection::~FeatureCollection() {
  int i;
  map< int, featureSet* >::iterator fiter;
  
  for (fiter = collection.begin(); fiter != collection.end(); fiter++) {
    closeFeatureSet(fiter->second);
  }
}

int FeatureCollection::addFeatureFile(const char* path, featureHeader* header, stegoContext* steg, featureSet* cleanSet) {
  int bin;
  featureSet *set;

  bin = (int) (header->prob*10000 + 0.5);
//   printf("using bin %i \n", bin);
  if (collection[bin] == 0) {
    if (cleanSet->header->slice_type == header->slice_type) { // maybe want more chekcs here
      set = openFeatureSet(path, steg);
      set->mask_vec = cleanSet->mask_vec;
      set->mask_vec = cleanSet->mask_vec;
      set->mask_counts = cleanSet->mask_counts;
      set->max_g    = cleanSet->max_g;
      set->min_g    = cleanSet->min_g;
      set->mu_g     = cleanSet->mu_g;
      set->mu_vec   = cleanSet->mu_vec;
      set->var_g    = cleanSet->var_g;
      
      set->prob = header->prob;
      set->id = bin;
      collection[bin] = set;
    }
  } else {
    if (collection[bin]->header->slice_type == header->slice_type)
      newFeatureFile(steg, collection[bin], path);
  }
}

int FeatureCollection::getNumSets() {
   return collection.size();
}

featureSet* FeatureCollection::getFeatureSet(int index) {
  return 0;
}

// int FeatureCollection::hasNext() {
//   if (current_set < num_sets) return 1;
//   else return 0;
// }

FeatureCollection::Iterator* FeatureCollection::iterator() {
  return new FeatureCollection::Iterator(this);
}

FeatureCollection::Iterator::Iterator(FeatureCollection* f) {
  fc = f;
  iter = f->collection.begin();
}

bool FeatureCollection::Iterator::hasNext() {
  return (iter != fc->collection.end());
}

featureSet* FeatureCollection::Iterator::next() {
  featureSet *set = iter->second;
  iter++;
  return set;
}


StegoModel::StegoModel() {
  int i, j;
//   features = new FeatureCollection();
//   current_view = 0;
  ranges = 0;
  steg = init_stego();
  cleanSet = 0;
  mc = 0;
  seenPaths = new set<string>();
//   for (i = 0; i < 10; i++) {
//     for (j = 0; j < 8; j++) {
//       collections[i][j] = 0;
//     }
//   }
}

StegoModel::~StegoModel() {
  if (mc != 0) closeMMD(*mc);
  close_stego(steg);
}

void StegoModel::estimateMus() {
  map< pair< int, int >, FeatureCollection* >::iterator fiter;
  FeatureCollection::Iterator *citer;
//   featureSet *stego;
  steg->features = cleanSet;
//   printf("activating something \n");
  startAction(steg->features);
//   printf("estimating some mu... \n");
  estimateMu(steg);
//   printf("done \n");
//   estimateSigma(steg);
//   qrHouseholder(steg);
  endAction(steg->features);
  for (fiter = collections.begin(); fiter != collections.end(); fiter++) {
    if (fiter->second != 0) {
//       printf("<%i, %i> \n", fiter->first.first, fiter->first.second);
      citer = fiter->second->iterator();
      while (citer->hasNext()) {
	steg->features = citer->next();
// 	printf("About to estimate mu with dim=%i, M=%i \n", steg->features->dim, steg->features->M);
	startAction(steg->features);
	estimateMu(steg);
	endAction(steg->features);
// 	progressChanged((double) i / (double) j);
      }
    }
  }
  steg->features = cleanSet;
  modelChanged();
}

// we don't want to give sets directly, rather indices inside the collection
double StegoModel::doMMD(featureSet *clean, featureSet *stego) {
//   mmdContext mc;
  if (mc == 0) {
    mc = (mmdContext*) malloc(sizeof(mmdContext));
    mc->clean = clean;
    mc->stego = stego;
    initMMD(steg, *mc);
    estimateGamma(steg, *mc);
  }
  mc->stego = stego;
//   mc->stego = mc->clean;
  estimateMMD(steg, *mc);
//   printf("used gamma: %g \n", mc->gamma);
  
  return mc->mmd;
}

void StegoModel::doMMDs() {
  map< pair< int, int >, FeatureCollection* >::iterator fiter;
  FeatureCollection::Iterator *citer;
//   featureSet *stego;
  
  mc = (mmdContext*) malloc(sizeof(mmdContext));
  mc->clean = cleanSet;
  startAction(mc->clean);
  initMMD(steg, *mc);
  estimateGamma(steg, *mc);
  
//   mc->stego = mc->clean;
//   estimateMMD(steg, *mc);
  
  for (fiter = collections.begin(); fiter != collections.end(); fiter++) {
    if (fiter->second != 0) {
      printf("<%i, %i> \n", fiter->first.first, fiter->first.second);
      citer = fiter->second->iterator();
      while (citer->hasNext()) {
        mc->stego = citer->next();
        printf("doing set %g \n", mc->stego->header->prob);
        startAction(mc->stego);
        estimateMMD(steg, *mc);
        mc->stego->mmd = mc->mmd;
        endAction(mc->stego);
      }
    }
  }
  endAction(mc->clean);
  closeMMD(*mc);
}

void StegoModel::runClassifier() {
  int dim = cleanSet->dim;
  std::cout << "dim = " << dim << std::endl;
  int read;
  double *vec_g, vec[dim];
  StegoClassifier sf(dim);
  
  CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
  startAction(cleanSet);
  for (int n = 0; n < cleanSet->M; n++) { // use while loop?
    readVectorRescaled(steg, cleanSet, vec_g);
    CUDA_CALL( cudaMemcpy(vec, vec_g, dim*sizeof(double), cudaMemcpyDefault));
    sf.addCleanVector(vec);
  }
  endAction(cleanSet);
}

void StegoModel::setFeatures(featureSet* set) {
  steg->features = set;
  modelChanged();
}


void StegoModel::addView(StegoView *view) {
  views.push_back(view);
//   current_view++;
}

void StegoModel::modelChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
//     printf("updateing some view \n");
    (*siter)->updateView();
  }
}

void StegoModel::collectionChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateCollection();
  }
}

void StegoModel::progressChanged(double p) {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateProgress(p);
  }
}

void StegoModel::openDirectory(const char* path) {
  int i, j, k;
  int num_sets = 0;
  int bin;
  char *str = (char*) malloc(512*sizeof(char));
  DIR *root = opendir(path);
  FILE *file;
  featureHeader header;
  struct dirent *entry;
  
  if (root == NULL) {
//     printf("root is NULL \n");
    return;
  }
//   printf("Root not NULL. \n");
  while ((entry = readdir(root)) != NULL) {
    if (strstr(entry->d_name, ".fv") != NULL) {
      num_sets++;
    }
  }
//   printf("about to open dir with %i feature files \n", num_sets);
  rewinddir(root);
  for(i = 0; i < num_sets; ) {
    entry = readdir(root);
    if (strstr(entry->d_name, ".fv") != NULL) {
      pathConcat(path, entry->d_name, str);
      openFile(str, i, num_sets, header);
//       file = fopen(str, "r");
//       readHeader(file, &header);
//       if (ranges == 0) {
// 	printf("first file is being added! \n");
// 	ranges = (int**) malloc(3*sizeof(int*));
// 	for (k = 0; k < 3; k++) {
// 	  ranges[k] = (int*) malloc(num_coefs[k]*sizeof(int));
// 	  for (j = 0; j < num_coefs[k]; j++)
// 	    ranges[k][j] = header.ranges[k][j];
// 	}
//       }
//       fclose(file);
// 
// //      printf("method: %i \n", header.method);
// //       printf("qp range: %i \n", header.qp_range);
//       if (header.method == 0) {
// 	if (cleanSet == 0) {
// 	  cleanSet = openFeatureSet(str, steg);
// 	} else {
// // 	  printf("adding new feature file to clean set \n");
// 	  newFeatureFile(steg, cleanSet, str);
// 	}
//       } else {
// 	if (collections[pair<int, int>(header.method,header.accept)] == 0) {
// 	  collections[pair<int, int>(header.method,header.accept)] = new FeatureCollection(&header);
// 	  printf("Created new collection for method %i and accept %i \n", header.method, header.accept);
// 	}
//         collections[pair<int, int>(header.method, header.accept)]->addFeatureFile(str, &header, steg, cleanSet);
//       }
//       i++;
//       progressChanged((double) i / (double) num_sets);
    }
  }
  collectionChanged();

  closedir(root);
  free(str);
}

int StegoModel::openFile(const char* path, int i, int num_sets, featureHeader &header) {
  int j, k;
  FILE *file;
//   featureHeader header;
  
  if (seenPaths->find(string(path)) != seenPaths->end()) {
//     printf("Das kenne ich doch schon %s :-@ \n", path);
    return 1;
  }
  
  file = fopen(path, "r");
  readHeader(file, &header);
  if (ranges == 0) {
    ranges = (int**) malloc(3*sizeof(int*));
    for (k = 0; k < 3; k++) {
      ranges[k] = (int*) malloc(num_coefs[k]*sizeof(int));
      for (j = 0; j < num_coefs[k]; j++)
        ranges[k][j] = header.ranges[k][j];
    }
    for (k = 0; k < 3; k++) {
      for (j = 0; j < num_coefs[k]; j++) {
	ranges[k][j] = 2 * (int) header.ranges[k][j] + 1;
// 	printf("ranges[%i][%i] = %i \n", k, j, ranges[k][j]);
      }
    }
  }
  fclose(file);

  if (header.method == 0) {
    if (cleanSet == 0) {
      cleanSet = openFeatureSet(path, steg);
    } else {
      newFeatureFile(steg, cleanSet, path);
    }
  } else {
    if (collections[pair<int, int>(header.method,header.accept)] == 0) {
      collections[pair<int, int>(header.method,header.accept)] = new FeatureCollection(&header);
    }
    collections[pair<int, int>(header.method, header.accept)]->addFeatureFile(path, &header, steg, cleanSet);
  }
  seenPaths->insert(string(path));
  progressChanged((double) i / (double) num_sets);
  
  return 0;
//   collectionChanged(); // maybe a bit inefficient to run this on every file, might be hundreds
}

StegoModel::Iterator::Iterator(StegoModel* m) {
  model = m;
  iter = m->collections.begin();
}

bool StegoModel::Iterator::hasNext() {
  return (iter != model->collections.end());
}

FeatureCollection* StegoModel::Iterator::next() {
  FeatureCollection *fc = iter->second;
  x = iter->first;
  iter++;
  return fc;
}

std::pair< int, int > StegoModel::Iterator::getX() {
  return x;
}

StegoModel::Iterator* StegoModel::iterator() {
  return new StegoModel::Iterator(this);
}

// FeatureCollection::Iterator* StegoModel::getFeatureIterator(int video_birate, int ppair, int method, int accept) {
// //   if (method < 0 || method >= 10) return 0;
// //   if (accept < 0 || accept >= 8)  return 0;
// //   if (collections[method][accept] != 0) {
// //     return collections[method][accept]->iterator();
// //   }
//   if (collections[pair<int, int>(method, accept)] != 0) {
//     return collections[pair<int, int>(method, accept)]->iterator();
//   }
//   return 0;
// }

featureSet* StegoModel::getCleanSet() {
  return cleanSet;
}

int StegoModel::getQPRange() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->header->qp_range;
}


int** StegoModel::getRanges() {
  return ranges;
}


int StegoModel::getDimension() {
  if (steg->features == NULL)
    return -1;
  return steg->features->dim;
}

int StegoModel::getHistDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->hist_dim;
}

int StegoModel::getPairDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->pair_dim;
}

int StegoModel::getUvsVDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->uvsv_dim;
}

double* StegoModel::getMaxVector() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->max_vec;
}

double* StegoModel::getMuVector() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->mu;
}

double* StegoModel::getQPHist() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->qp_vec;
}

// FeatureCollection* StegoModel::getCollection() {
// //   return fcol; // fcol
//   return 0;
// }

int StegoModel::getSigmaDim() {
  if (mc == 0) return -1;
  return mc->cache;
}

double* StegoModel::getSigma() {
  if (steg->features == 0) 
    return 0;
  return steg->features->gauss->sigma;
//   if (mc == 0) return 0;
//   return mc->results;
}

double* StegoModel::getDiag() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->qr_diag;
}

// StegoModel::Iterator* StegoModel::getIterator() {
//   return new Iterator(this);
// }
// 
// StegoModel::Iterator::Iterator(StegoModel* model) {
//   level = 0;
//   level1 = -1;
//   level3 = -1;
// //   level0iter = model->collections.iterator();
//   next();
// }
// 
// int StegoModel::Iterator::getLevel0() {
// //   return *level0iter->first;
// }
// 
// int StegoModel::Iterator::getLevel1() {
//   return level1;
// }
// 
// int StegoModel::Iterator::getLevel2() {
// //   return level2iter->first;
// }
// 
// int StegoModel::Iterator::getLevel3() {
//   return level3;
// }
// 
// bool StegoModel::Iterator::hasNext() {
//   return true;
// }
// 
// void StegoModel::Iterator::next() {
//   int i = 3;
//   
// }

// // FeatureCollection::Iterator* StegoModel::Iterator::nextCollectionIterator() {
// //   FeatureCollection::Iterator* iter = level3iter->iterator();
// //   next();
// //   return iter;
// // }

