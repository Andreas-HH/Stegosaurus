#include "Stegosaurus.h"


void initMMD(stegoContext *steg, mmdContext& mc) {
  int i;
  int dim = mc.clean->dim;
  
  mc.n = mc.clean->M;
  mc.cache   = min(mc.n, (int) (steg->gpu_c->doublesOnGPU / (long) dim / 4l)); // 4 input buffers
  
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_right_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_right_g, dim*mc.cache*sizeof(double)));
//   CUDA_CALL( cudaHostAlloc(&mc.vectors, dim*mc.cache*sizeof(double), cudaHostAllocDefault));
//   CUDA_CALL( cudaHostAlloc(&mc.stego_vectors, dim*mc.cache*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mc.v_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&mc.v2_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.temp_g, dim*sizeof(double)));  
  
//   for (i = 0; i < mc.n; i++) {
//     readVectorL2(steg, mc.clean, mc.clean_vectors+i*dim);
//     readVectorL2(steg, mc.stego, mc.stego_vectors+i*dim);
//   }
//   stegoRewind(mc.clean);
//   stegoRewind(mc.stego);
}

void closeMMD(mmdContext& mc) {
  CUDA_CALL( cudaFree(mc.clean_vectors_down_g));
  CUDA_CALL( cudaFree(mc.clean_vectors_right_g));
  CUDA_CALL( cudaFree(mc.stego_vectors_down_g));
  CUDA_CALL( cudaFree(mc.stego_vectors_right_g));
  CUDA_CALL( cudaFree(mc.v_g));
  CUDA_CALL( cudaFree(mc.temp_g));
//   CUDA_CALL( cudaFreeHost(mc.vectors));
//   CUDA_CALL( cudaFree(mc.vectors_g));
//   CUDA_CALL( cudaFree(mc.v_g));
//   CUDA_CALL( cudaFree(mc.temp_g));
}

void estimateGamma(stegoContext *steg, mmdContext& mc) {
  int i, j;
  int blockwidth, pos;
//   int iterations;
  double norm;
//   double *vectors = mc->vectors;
//   double *v1 = mc->v_g;//, *v2 = mc->v2_g;
  double min1 = -1;
  featureSet *cleanSet = mc.clean;
//   int r;
  int M = mc.n;
  int dim = cleanSet->dim;
  priority_queue< double > q;
  
//   printf("cache size: %i, M= %i \n", mc.cache, M);
  // could use 4 times the buffer size here
//   iterations = (M+mc.cache-1)/mc.cache;
//   printf("will need %i iterations, caching %i vectors \n", iterations, mc.cache);
  for (pos = 0; pos < M; pos += mc.cache) {
//     printf("nw in position: %i / %i \n", pos, M);
    blockwidth = min(mc.cache, M-pos);
    for (i = 0; i < blockwidth; i++) {
      readVectorRescaled(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
    }
    for (i = 0; i < M; i++) {
//       if (i % 100 == 0) printf("i = %i, blockwidth = %i \n", i, blockwidth);
      readVectorRescaled(steg, mc.stego, mc.v_g);
      for (j = 0; j < blockwidth; j++) {
	if (pos+j == i) continue;
	CUDA_CALL( cudaMemcpy(mc.temp_g, mc.v_g, dim*sizeof(double), cudaMemcpyDeviceToDevice));
	// maybe I should write my own subtraction kernel to abvoid the unnecessary multiplication by 1
	CUBLAS_CALL( cublasDaxpy(steg->gpu_c->handle, dim, &min1, mc.clean_vectors_down_g+j*dim, 1, mc.temp_g, 1));
        CUBLAS_CALL( cublasDdot(steg->gpu_c->handle, dim, mc.temp_g, 1, mc.temp_g, 1, &norm));
        q.push(norm);
      }
    }
    stegoRewind(mc.stego);
  }
  stegoRewind(mc.clean);
  
  printf("queue size: %i, M = %i, expcted size: %i \n", q.size(), M, M*(M-1));
  for (i = 0; i < M*(M-1)/2; i++) {
    q.pop();
  }
  printf("median: %g => gamma = %g , queue size: %i \n", q.top(), 1./q.top(), q.size());
  mc.gamma = 1./q.top();
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

void estimateMMD(stegoContext *steg, mmdContext& mc) {
  int i, j, k;
  int bw_x, bw_y;
  long pos_x, pos_y;
//   int cache = mc->cache;
//   int blockwidth_g, pos_g, cache_g = mc->cache_g;
//   int lessThanOne = 0, greaterThanOne = 0;
//   double *vectors;
  int dim = mc.clean->dim;
  int M = mc.n;
  double mmd = 0.;
  double temp;
  
  for (pos_x = 0l; pos_x < (long) M; pos_x += mc.cache) {
    bw_x = min(mc.cache, M-(int)pos_x);
    printf("clean [%i] \n", pos_x);
    jumpToVector(mc.clean, pos_x);
    printf("stego [%i] \n", pos_x);
    jumpToVector(mc.stego, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
    } // better to one after the other, might be quicker
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.stego, mc.stego_vectors_down_g + i*dim);
    }
    for (pos_y = 0l; pos_y < (long) (M+mc.cache-1); pos_y += mc.cache) {
      bw_y = min(mc.cache, M-(int)pos_y);
      printf(" clean [%i] \n", pos_y);
      jumpToVector(mc.clean, pos_y);
      printf(" stego [%i] \n", pos_y);
      jumpToVector(mc.stego, pos_y);
      for (i = 0; i < bw_x; i++) {
        readVectorRescaled(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
      }
      for (i = 0; i < bw_x; i++) {
        readVectorRescaled(steg, mc.stego, mc.stego_vectors_right_g + i*dim);
      }
      for (j = 0; j < bw_x; j++) {
	for (k = 0; k < bw_y; k++) {
	  if (pos_x+j == pos_y + k)  continue;
	  mmd += applyKernel(steg, mc.gamma, dim, mc.clean_vectors_down_g+j*dim, mc.clean_vectors_right_g+k*dim, mc.temp_g);
	  mmd += applyKernel(steg, mc.gamma, dim, mc.stego_vectors_down_g+j*dim, mc.stego_vectors_right_g+k*dim, mc.temp_g);
	  mmd -= 2.*applyKernel(steg, mc.gamma, dim, mc.clean_vectors_down_g+j*dim, mc.stego_vectors_right_g+k*dim, mc.temp_g);
	}
      }
    }
//     stegoRewind(mc.stego);
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
  printf("have some mmd: %g \n", mmd);
  mc.mmd = 1./((double) M*(M-1)) * mmd;
}