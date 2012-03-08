#include "Stegosaurus.h"


void loadVectorsMMD(stegoContext *steg, mmdContext& mc) {
  int i;
  int dim = mc.clean->dim;
  
  mc.n = mc.clean->M;
  
  CUDA_CALL( cudaHostAlloc(&mc.clean_vectors, dim*mc.n*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&mc.stego_vectors, dim*mc.n*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mc.v1_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.v2_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.temp_g, dim*sizeof(double)));  
  
  for (i = 0; i < mc.n; i++) {
    readVectorL2(steg, mc.clean, mc.clean_vectors+i*dim);
    readVectorL2(steg, mc.stego, mc.stego_vectors+i*dim);
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
}

void estimateGamma(stegoContext *steg, mmdContext *mc) {
  int i, j;
  double norm;
  double *vectors = mc->clean_vectors;
  double *v1 = mc->v1_g, *v2 = mc->v2_g;
  double min1 = -1;
  featureSet *cleanSet = mc->clean;
//   int r;
  int M = mc->n;
  int dim = cleanSet->dim;
  priority_queue< double > q;
  
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (i == j) continue;
      CUDA_CALL( cudaMemcpy(v1, vectors+i*dim, dim*sizeof(double), cudaMemcpyDeviceToDevice));
      CUBLAS_CALL( cublasDaxpy(steg->gpu_c->handle, dim, &min1, vectors+j*dim, 1, v1, 1));
      CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, dim, v1, 1, v1, 1, &norm));
      q.push(norm);
    }
  }
  printf("queue size: %i, M = %i, expcted size: %i \n", q.size(), M, M*(M-1));
  for (i = 0; i < M*(M-1)/2; i++) {
    q.pop();
  }
  printf("median: %g => gamma = %g , queue size: %i \n", q.top(), 1./q.top(), q.size());
  mc->gamma = 1./q.top();
}

// probably want something block- rather than tuple-wise later!
double applyKernel(stegoContext *steg, double gamma, int dim, double *v1_g, double *v2_g, double *temp_g) {
  double norm; // use dotp instead?
  double min1 = -1;
  
  CUDA_CALL( cudaMemcpy(temp_g, v2_g, dim*sizeof(double), cudaMemcpyDeviceToDevice));
  CUBLAS_CALL( cublasDaxpy(steg->gpu_c->handle, dim, &min1, v1_g, 1, temp_g, 1));
  CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, dim, temp_g, 1, temp_g, 1, &norm));
  
  return exp(-1*gamma*norm);
}

void estimateMMD(stegoContext *steg, mmdContext *mc) {
  int i, j, k;
  int blockwidth, pos;
  int lessThanOne = 0, greaterThanOne = 0;
  double *vectors;
  int dim = mc->clean->dim;
  int max_on_gpu = 10000;
  int M = mc->n;
  double mmd = 0.;
  double temp;
  
  CUDA_CALL( cudaHostAlloc(&vectors, mc->clean->dim*M*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mc->vectors_g, M*dim*sizeof(double)));
  
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (i == j) continue;
      temp = applyKernel(steg, mc->gamma, dim, mc->clean_vectors+i*dim, mc->clean_vectors+j*dim, mc->temp_g);
      mmd += temp;
      if (temp < exp(-1.)) greaterThanOne++;
      else lessThanOne++;
    }
  }
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (i == j) continue;
      temp = applyKernel(steg, mc->gamma, dim, mc->stego_vectors+i*dim, mc->stego_vectors+j*dim, mc->temp_g);
      mmd += temp;
      if (temp < exp(-1.)) greaterThanOne++;
      else lessThanOne++;
    }
  }
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (i == j) continue;
      temp = applyKernel(steg, mc->gamma, dim, mc->clean_vectors+i*dim, mc->stego_vectors+j*dim, mc->temp_g);
      mmd -= 2*temp;
      if (temp < exp(-1.)) greaterThanOne++;
      else lessThanOne++;
    }
  }
  printf("current mmd: %g, %i <= 1, %i > 1 \n", mmd, lessThanOne, greaterThanOne);
  mc->mmd = mmd;
}