#include "Stegosaurus.h"


__global__ void compareMax(int dim, double *current_max, double *new_features) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    if (current_max[idx] < new_features[idx])
      current_max[idx] = new_features[idx];
  }
}

__global__ void subtract(int dim, double *vec, double *mu) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    vec[idx] -= mu[idx];
  }
}

/*int computeQPHistogram(stegoContext *steg, double *mu_g, int qp_range, double *result) {
  int i;
  double alpha = 100./((double) steg->features->dim);
  
  for (i = 0; i < qp_range; i++) {
    cublasDasum_v2(steg->gpu_c->handle, steg->features->dim/qp_range, steg->features->mu_g+i*steg->features->dim/qp_range, 1, result+i);
  }
  cublasDscal_v2(steg->gpu_c->handle, 20, &alpha, result, 1);
  
  return 0;
}*/

int estimateMu(stegoContext *steg) {
  int i;
  int read;
//   double sum;
  featureSet *fs = steg->features;
  int dim = fs->dim;
  int tpb = steg->gpu_c->threads_per_block;
//   double zero = 0.;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
  double *current_feature;
  double *mu_g, *max_g, *vec_g;
  
  if (fs->gauss->mu != NULL) return -1;
  
  CUDA_CALL( cudaHostAlloc(&current_feature, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mu_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&qp_g, 20*sizeof(double)));
  CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&max_g, dim*sizeof(double)));
  
  fs->gauss->mu = (double*) malloc(dim*sizeof(double));
  for (i = 0; i < dim; i++)
    current_feature[i] = 0.;

  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, mu_g, 1));
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, max_g, 1));

//   printf("uploaded stuff \n");
  for (i = 0; i < fs->M; i++) {
     read = readVector(fs, current_feature);
     if (read != dim) printf("read something wrong: %i \n", read);
     CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
     cublasDaxpy(*handle, dim, &(fs->divM), vec_g, 1, mu_g, 1);     // put divM on gpu?
     compareMax<<<BLOCKS(dim,tpb),tpb>>>(dim, max_g, vec_g);
  }
//   computeQPHistogram(steg, fs->mu_g, 20, fs->qp_g);
//   CUBLAS_CALL( cublasGetVector(20, sizeof(double), fs->qp_g, 1, fs->qp_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), mu_g, 1, fs->gauss->mu, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), max_g, 1, fs->max_vec, 1));
  stegoRewind(fs);
  
  for (i = 0; i < QP_RANGE; i++) {
    fs->qp_vec[i] = 0;
    for (int j = 0; j < dim/QP_RANGE; j++) {
      fs->qp_vec[i] += fs->gauss->mu[i*dim/QP_RANGE+j];
    }
    fs->qp_vec[i] *= 100./((double) dim);
  }
  
  CUDA_CALL( cudaFreeHost(current_feature));
  CUDA_CALL( cudaFree(mu_g));
//   CUDA_CALL( cudaFree(qp_g));
  CUDA_CALL( cudaFree(vec_g));
  CUDA_CALL( cudaFree(max_g));
  
  return 0;
}

int estimateSigma(stegoContext *steg) {
  int i, j;
  int read;
  int blockwidth, posInSigma;   // size of current block on gpu
  featureSet *fs = steg->features;
  myGaussian *gauss = fs->gauss;
//   double minus_one = -1;
  int dim = fs->dim;
  int tpb = steg->gpu_c->threads_per_block;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
  double *current_feature;
  double *mu_g, *vec_g, *sigma_g;
  
  if (!(gauss->mu != NULL && gauss->sigma == NULL)) return -1;
  
//   gauss->sigma = (double*) malloc(dim*dim*sizeof(double));
  CUDA_CALL( cudaHostAlloc(&(gauss->sigma), dim*dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&current_feature, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mu_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&sigma_g, dim*fs->gpu_matrix_width*sizeof(double)));
  
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), gauss->mu, 1, mu_g, 1));
//   CUBLAS_CALL( cublasDscal(*handle, dim, &(minus_one), mu_g, 1));
  
  for (i = 0; i < dim*dim; i++) {
    gauss->sigma[i] = 0;
  }
  
  for (j = 0; j < (dim+fs->gpu_matrix_width-1)/fs->gpu_matrix_width; j++) { // (dim+fs->gpu_matrix_width-1)/fs->gpu_matrix_width
    posInSigma = j*fs->gpu_matrix_width;
    blockwidth = min(fs->gpu_matrix_width, dim-posInSigma);
    CUBLAS_CALL( cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->sigma+posInSigma*dim, dim, sigma_g, dim)); // add on cpu??
    printf("about to calc sigma, blockwidth=%i, dim=%i, pos=%i mw=%i\n", blockwidth, dim, posInSigma, fs->gpu_matrix_width);
    for (i = 0; i < fs->M; i++) { // fs->M
      read = readVector(fs, current_feature);
      if (read != dim) printf("read something wrong: %i \n", read);
      CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
      subtract<<<BLOCKS(dim,tpb),tpb>>>(dim, vec_g, mu_g);
      
//       cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->sigma+j*fs->gpu_matrix_width*dim, dim, sigma_g, dim); // add on cpu??
      
      CUBLAS_CALL( cublasDger(*handle, dim, blockwidth, &(fs->divM), vec_g, 1, vec_g+j*fs->gpu_matrix_width, 1, sigma_g, dim));
      if (i%100 == 0) printf("%i \n", i);
//       cublasGetMatrix(dim, blockwidth, sizeof(double), sigma_g, dim, gauss->sigma+j*fs->gpu_matrix_width*dim, dim); // only do this with multiple blocks
    }
    stegoRewind(fs);
    CUBLAS_CALL( cublasGetMatrix(dim, blockwidth, sizeof(double), sigma_g, dim, gauss->sigma+posInSigma*dim, dim)); // only do this with multiple blocks
  }
  
  CUDA_CALL( cudaFree(sigma_g));
  CUDA_CALL( cudaFree(mu_g));
  CUDA_CALL( cudaFree(vec_g));
  CUDA_CALL( cudaFreeHost(current_feature));
  
  return 0;
}

