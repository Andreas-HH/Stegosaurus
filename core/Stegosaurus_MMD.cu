#include "Stegosaurus.h"

void estimateGamma(stegoContext *steg, mmdContext *mc) {
  int i;
  double norm;
  double *vectors;
  double *v1, *v2;
  double min1 = -1;
  featureSet *stegoSet = mc->stego;
  int r;
  int M = min(stegoSet->M, mc->clean->M);
  int dim = stegoSet->dim;
  
//   steg->features = stegoSet;
  mc->n = M;
  CUDA_CALL( cudaHostAlloc(&vectors, stegoSet->dim*stegoSet->M*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&v1, stegoSet->dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&v2, stegoSet->dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc->temp_g, stegoSet->dim*sizeof(double)));
  
  for (i = 0; i < M; i++) {
    readVector(steg, stegoSet, vectors+i*dim);
  }
  stegoRewind(stegoSet);
  mc->stego_vectors = vectors;
  mc->v1_g = v1;
  mc->v2_g = v2;
  
  for (i = 0; i < 10; i++) {
    r = rand() % M;
    CUBLAS_CALL( cublasSetVector(stegoSet->dim, sizeof(double), vectors + r*dim, 1, v1, 1));
    CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, dim, v1, 1, &norm));
    norm = 1./norm;
    CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, dim, &norm, v1, 1));
    r = rand() % M;
    CUBLAS_CALL( cublasSetVector(stegoSet->dim, sizeof(double), vectors + r*dim, 1, v2, 1));
    CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, dim, v2, 1, &norm));
    norm = 1./norm;
    CUBLAS_CALL( cublasDscal(steg->gpu_c->handle, dim, &norm, v2, 1));
    printf("applying kernel to something: %g ", applyKernel(steg, 1.5, dim, v1, v2, mc->temp_g));
    CUBLAS_CALL( cublasDaxpy(steg->gpu_c->handle, dim, &min1, v2, 1, v1, 1));
    CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, dim, v1, 1, &norm));
    printf("r=%i, norm: %g \n", r, norm);
  //   r = (rand()*stegoSet->M)/RAND_MAX;
  //   CUBLAS_CALL( cublasSetVector(stegoSet->dim, sizeof(double), vectors + r*stegoSet->dim*sizeof(double), 1, v1, 1));
  }
  mc->gamma = 1.5;  // needs to be done properly
  
//   CUDA_CALL( cudaFree(v1));
//   CUDA_CALL( cudaFree(v2));
//   CUDA_CALL( cudaFreeHost(vectors));
}

double applyKernel(stegoContext *steg, double gamma, int dim, double *v1_g, double *v2_g, double *temp_g) {
  double norm; // use dotp instead?
  double min1 = -1;
  
  CUDA_CALL( cudaMemcpy(temp_g, v2_g, dim*sizeof(double), cudaMemcpyDeviceToDevice));
  CUBLAS_CALL( cublasDaxpy(steg->gpu_c->handle, dim, &min1, v1_g, 1, temp_g, 1));
  CUBLAS_CALL( cublasDnrm2(steg->gpu_c->handle, dim, temp_g, 1, &norm));
  
  printf("norm: %g \n", norm);
  
  return exp(-1*gamma*norm*norm);
}

void estimateMMD(stegoContext *steg, mmdContext *mc) {
  int i, j, k;
  int blockwidth, pos;
  double *vectors;
  int dim = mc->clean->dim;
  int max_on_gpu = 10000;
  int M = mc->n;
  double mmd = 0.;
  
  CUDA_CALL( cudaHostAlloc(&vectors, mc->clean->dim*M*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mc->vectors_g, M*dim*sizeof(double)));
  
  for (i = 0; i < M; i++) {
    readVector(steg, mc->clean, vectors+i*dim);
  }
  stegoRewind(mc->clean);
  
  for (i = 0; i < (M+max_on_gpu-1)/max_on_gpu; i++) {
    pos = i * M;
    blockwidth = min(M-pos, max_on_gpu);
//     printf("bw: %i \n", blockwidth);
    // upload a block of x vectors
    CUBLAS_CALL( cublasSetVector(blockwidth*dim, sizeof(double), vectors+pos*dim, 1, mc->vectors_g, 1));
    for (j = 0; j < M; j++) {
      // upload a y vector
      CUBLAS_CALL( cublasSetVector(dim, sizeof(double), mc->stego_vectors+j*dim, 1, mc->v2_g, 1));
      for (k = 0; k < blockwidth; k++) {
	if (pos + k != j) {
	  // apply kernel on x_k and y_j
	  mmd = applyKernel(steg, mc->gamma, dim, mc->vectors_g+k*dim, mc->v2_g, mc->temp_g); // vectors_g+k*dim
	  printf("applying x_%i to y_%i with gamma=%g: %g \n", pos+k, j, mc->gamma, mmd);
// 	  mmd -= 2*applyKernel(steg, mc->gamma, dim, mc->vectors_g+k*dim, mc->v2_g, mc->temp_g);
	}
      }
    }
  }
  
  printf("negative part: %g \n", mmd);
}