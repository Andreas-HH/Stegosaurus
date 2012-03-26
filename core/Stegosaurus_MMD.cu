#include "Stegosaurus.h"


__global__ void gammaKernel(int dim, int cache, int offset, int steps, int bw_x, int bw_y, double *down_g, double *right_g, double *results) {
  int i;
  int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
  int idx_y = threadIdx.y + blockIdx.y*blockDim.y;
  int dx = idx_x*dim + offset;
  int dy = idx_y*dim + offset;
  double temp;
  double current_sum = 0.;
  
  if (idx_x < bw_x && idx_y < bw_y) {
    for (i = 0; i < steps; i++) {
      temp = down_g[dy + i] - right_g[dx + i];
      current_sum += temp * temp;
    }
//     results[idx_y + bw_y*idx_x] += temp * temp;//current_sums[idx_s];
    results[idx_y + cache*idx_x] += current_sum;
  }
}

__global__ void mmdKernel(double minus_gamma, double *cvc_g, double *cvs_g, double *svs_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
//   double cvc = minus_gamma * cvc_g[idx] + 1.;
  double cvc = exp(minus_gamma * cvc_g[idx]);
  double cvs = exp(minus_gamma * cvs_g[idx]);
  double svs = exp(minus_gamma * svs_g[idx]);
  
//   if (cvc < 0) cvc = -1.;
//   if (cvc > 0) cvc = 1.;
  
  cvc_g[idx] = cvc + svs - 2*cvs;
}

void initMMD(stegoContext *steg, mmdContext& mc) {
  int dim = mc.clean->dim;
//   int tpb = steg->gpu_c->threads_per_block;
  
  mc.n = mc.clean->M;
  mc.kernel_blockwidth = (int) sqrt(steg->gpu_c->threads_per_block);
  mc.cache   = min(mc.n, (int) (sqrt(steg->gpu_c->doublesOnGPU / 3l + (long) SQUARE(dim) * 4l / 9l) - (long)dim * 2l / 3l));
  mc.kernel_gridwidth = (mc.cache + mc.kernel_blockwidth-1)/mc.kernel_blockwidth;
//   printf("cache: %i, kbw: %i \n", mc.cache, mc.kernel_blockwidth);
  
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_right_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_right_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_c_vs_c_g, mc.cache*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_c_vs_s_g, mc.cache*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_s_vs_s_g, mc.cache*mc.cache*sizeof(double)));
  CUDA_CALL( cudaHostAlloc(&mc.results, mc.cache*mc.cache*sizeof(double), cudaHostAllocDefault));
}

void closeMMD(mmdContext& mc) {
  CUDA_CALL( cudaFree(mc.clean_vectors_down_g));
  CUDA_CALL( cudaFree(mc.clean_vectors_right_g));
  CUDA_CALL( cudaFree(mc.stego_vectors_down_g));
  CUDA_CALL( cudaFree(mc.stego_vectors_right_g));
  CUDA_CALL( cudaFree(mc.results_c_vs_c_g));
  CUDA_CALL( cudaFree(mc.results_c_vs_s_g));
  CUDA_CALL( cudaFree(mc.results_s_vs_s_g));
  CUDA_CALL( cudaFreeHost(mc.results));
}

void launchGammaKernel(mmdContext& mc, int dim, int bw_x, int bw_y, double* down_g, double* right_g, double* results_g) {
  int i;
  int step_size = 128;
  dim3 grid, block;
  
  grid = dim3(BLOCKS(bw_x, mc.kernel_blockwidth), BLOCKS(bw_y, mc.kernel_blockwidth));    
  block = dim3(mc.kernel_blockwidth, mc.kernel_blockwidth);
  for (i = 0; i < dim; i += step_size) {
//     if (i % 100 == 0) printf("i = %i / %i \n", i, dim);
    gammaKernel<<<grid,block>>>(dim, mc.cache, i, min(step_size, dim - i), bw_x, bw_y, down_g, right_g, results_g);
    cudaThreadSynchronize();
  }
}

void estimateGamma(stegoContext *steg, mmdContext& mc) {
  int i, j;
  int bw_x, bw_y, pos_x, pos_y;
  int tpb = steg->gpu_c->threads_per_block;
  featureSet *cleanSet = mc.clean;
  long M = (long) mc.n;
  long l;
  int dim = cleanSet->dim;
  priority_queue< double > q;
//   time_t start;
  
  for (pos_x = 0l; pos_x < M; pos_x += mc.cache) {
    bw_x = min(mc.cache, (int) (M-pos_x));
    jumpToVector(mc.clean, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorL2(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
    }
    for (pos_y = pos_x; pos_y < M; pos_y += mc.cache) {
      bw_y = min(mc.cache, (int) (M-pos_y));
      jumpToVector(mc.clean, pos_y);
      for (i = 0; i < bw_y; i++) {
        readVectorL2(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
      }
      printf("launching kernel with parameters (%i, %i), (%i, %i), bw_x = %i, bw_y = %i", BLOCKS(bw_x, mc.kernel_blockwidth), BLOCKS(bw_y, mc.kernel_blockwidth), mc.kernel_blockwidth, mc.kernel_blockwidth, bw_x, bw_y);
      initDArray(mc.results_c_vs_c_g, SQUARE(mc.cache), tpb, 0.);
//       start = time(NULL);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.clean_vectors_right_g, mc.results_c_vs_c_g);
//       printf(" ... took %is \n", time(NULL)-start);
      CUDA_CALL( cudaMemcpy(mc.results, mc.results_c_vs_c_g, SQUARE(mc.cache)*sizeof(double), cudaMemcpyDeviceToHost));
//       CUBLAS_CALL( cublasGetVector(SQUARE(mc.cache), sizeof(double), mc.results_c_vs_c_g, 1, mc.results, 1));
      for (i = 0; i < bw_x; i++) {
	for (j = 0; j < bw_y; j++) {
// 	  if (pos_x + i == pos_y + j) continue;
          if (pos_x + i < pos_y + j)  {
	    q.push(mc.results[j + i*mc.cache]);
	  }
	}
      }
    }
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
  
  printf("queue size: %i, M = %i, expcted size: %i \n", q.size(), M, M*(M-1)/2);
  for (l = 0; l < M*(M-1)/4l; l++) {
    q.pop();
  }
  printf("median: %g => gamma = %g , queue size: %i \n", q.top(), 1./q.top(), q.size());
  mc.gamma = 1./q.top();
}

void estimateMMD(stegoContext *steg, mmdContext& mc) {
  int i, j, k;
  int bw_x, bw_y;
  long pos_x, pos_y;
  int tpb = steg->gpu_c->threads_per_block;
  int dim = mc.clean->dim;
  int M = mc.n;
  double mmd = 0.;
  time_t start = time(NULL);
  
  for (pos_x = 0l; pos_x < (long) M; pos_x += mc.cache) {
    bw_x = min(mc.cache, M-(int)pos_x);
    jumpToVector(mc.clean, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorL2(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
    }
    jumpToVector(mc.stego, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorL2(steg, mc.stego, mc.stego_vectors_right_g + i*dim);
    }
    for (pos_y = 0l; pos_y < (long) M; pos_y += mc.cache) {
      bw_y = min(mc.cache, M-(int)pos_y);
      jumpToVector(mc.clean, pos_y);
      for (i = 0; i < bw_y; i++) {
        readVectorL2(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
      }
      jumpToVector(mc.stego, pos_y);
      for (i = 0; i < bw_y; i++) {
        readVectorL2(steg, mc.stego, mc.stego_vectors_down_g + i*dim);
      }
//       printf("launching kernel with parameters (%i, %i), (%i, %i), bw_x = %i, bw_y = %i", BLOCKS(bw_x, mc.kernel_blockwidth), BLOCKS(bw_y, mc.kernel_blockwidth), mc.kernel_blockwidth, mc.kernel_blockwidth, bw_x, bw_y);
      initDArray(mc.results_c_vs_c_g, SQUARE(mc.cache), tpb, 0.);
      initDArray(mc.results_c_vs_s_g, SQUARE(mc.cache), tpb, 0.);
      initDArray(mc.results_s_vs_s_g, SQUARE(mc.cache), tpb, 0.);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.clean_vectors_right_g, mc.results_c_vs_c_g);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.stego_vectors_right_g, mc.results_c_vs_s_g);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.stego_vectors_down_g, mc.stego_vectors_right_g, mc.results_s_vs_s_g);
      mmdKernel<<<BLOCKS( mc.cache*mc.cache, tpb), tpb>>>(-1.*mc.gamma, mc.results_c_vs_c_g, mc.results_c_vs_s_g, mc.results_s_vs_s_g);
//       cudaThreadSynchronize();
      CUBLAS_CALL( cublasGetVector(SQUARE(mc.cache), sizeof(double), mc.results_c_vs_c_g, 1, mc.results, 1));
//       CUBLAS_CALL( cublasSetVector(SQUARE(mc.cache), sizeof(double), mc.results, 1, mc.results_c_vs_c_g, 1));
      for (i = 0; i < bw_x; i++) {
	for (j = 0; j < bw_y; j++) {
	  if (pos_x + i == pos_y + j) continue;
          mmd += mc.results[j + i*mc.cache];
	}
      }
    }
//     break;
//     stegoRewind(mc.stego);
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
  mmd /= (double) (M * (M-1));
  printf("have some mmd: %f [%is]\n", mmd, time(NULL)-start);
  mc.mmd = mmd;
}