#include "Stegosaurus.h"


__global__ void gammaKernel(int dim, uint64_t cache, int offset, int steps, uint64_t bw_x, uint64_t bw_y, double *down_g, double *right_g, double *results) {
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
    results[idx_y + cache*idx_x] += current_sum;
  }
}

__global__ void mmdKernel(double minus_gamma, double *cvc_g, double *cvs_g, double *svs_g) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  double cvc = exp(minus_gamma * cvc_g[idx]);
  double cvs = exp(minus_gamma * cvs_g[idx]);
  double svs = exp(minus_gamma * svs_g[idx]);
  
  cvc_g[idx] = cvc + svs - 2*cvs;
}

void initMMD(stegoContext *steg, mmdContext& mc) {
  int dim = mc.clean->dim;
  
  mc.n = mc.clean->M;
  mc.kernel_blockwidth = (int) sqrt(steg->gpu_c->threads_per_block);
  mc.cache   = min(mc.n, (uint64_t) (sqrt(steg->gpu_c->doublesOnGPU / 3l + (long) SQUARE(dim) * 4l / 9l) - (long)dim * 2l / 3l));
  mc.kernel_gridwidth = (mc.cache + mc.kernel_blockwidth-1)/mc.kernel_blockwidth;
  printf("cache: %d, kbw: %i \n", mc.cache, mc.kernel_blockwidth);
  
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

void launchGammaKernel(mmdContext& mc, int dim, uint64_t bw_x, uint64_t bw_y, double* down_g, double* right_g, double* results_g) {
  int i;
  int step_size = 128;
  dim3 grid, block;
  
  grid = dim3(BLOCKS(bw_x, mc.kernel_blockwidth), BLOCKS(bw_y, mc.kernel_blockwidth));    
  block = dim3(mc.kernel_blockwidth, mc.kernel_blockwidth);
  for (i = 0; i < dim; i += step_size) {
    gammaKernel<<<grid,block>>>(dim, mc.cache, i, min(step_size, dim - i), bw_x, bw_y, down_g, right_g, results_g);
    cudaThreadSynchronize();
  }
}

void estimateGamma(stegoContext *steg, mmdContext& mc) {
  uint64_t i, j;
  uint64_t bw_x, bw_y, pos_x, pos_y;
  int tpb = steg->gpu_c->threads_per_block;
  featureSet *cleanSet = mc.clean;
  uint64_t M = mc.n;
  uint64_t l;
  int dim = cleanSet->dim;
  priority_queue< double > q;
  
  for (pos_x = 0ull; pos_x < M; pos_x += mc.cache) {
    bw_x = min(mc.cache, M-pos_x);
    jumpToVector(mc.clean, pos_x);
    for (i = 0ull; i < bw_x; i++) {
      readVectorRescaled(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
    }
    for (pos_y = pos_x; pos_y < M; pos_y += mc.cache) {
      bw_y = min(mc.cache, M-pos_y);
      jumpToVector(mc.clean, pos_y);
      for (i = 0ull; i < bw_y; i++) {
        readVectorRescaled(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
      }
      initDArray(mc.results_c_vs_c_g, SQUARE(mc.cache), tpb, 0.);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.clean_vectors_right_g, mc.results_c_vs_c_g);
      CUDA_CALL( cudaMemcpy(mc.results, mc.results_c_vs_c_g, SQUARE(mc.cache)*sizeof(double), cudaMemcpyDeviceToHost));
      cudaThreadSynchronize();
      for (i = 0ull; i < bw_x; i++) {
        for (j = 0ull; j < bw_y; j++) {
              if (pos_x + i < pos_y + j)  {
            q.push(mc.results[j + i*mc.cache]);
          }
        }
      }
    }
  }
  stegoRewind(mc.clean);
  
  printf("queue size: %i, M = %i, expected size: %i \n", q.size(), M, M*(M-1ull)/2ull);
  for (l = 0; l < M*(M-1ull)/4ull; l++) {
    q.pop();
  }
  printf("median: %g => gamma = %g , queue size: %i \n", q.top(), 1./q.top(), q.size());
  mc.gamma = 1./q.top();
}

void estimateMMD(stegoContext *steg, mmdContext& mc) {
  uint64_t i, j, k;
  uint64_t bw_x, bw_y;
  uint64_t pos_x, pos_y;
  int tpb = steg->gpu_c->threads_per_block;
  int dim = mc.clean->dim;
  uint64_t M = mc.n;
  double mmd = 0.;
  time_t start = time(NULL);
  
  for (pos_x = 0ull; pos_x < M; pos_x += mc.cache) {
    bw_x = min(mc.cache, M-pos_x);
    jumpToVector(mc.clean, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
    }
    jumpToVector(mc.stego, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.stego, mc.stego_vectors_right_g + i*dim);
    }
    for (pos_y = 0ull; pos_y < M; pos_y += mc.cache) {
      bw_y = min(mc.cache, M-pos_y);
      jumpToVector(mc.clean, pos_y);
      for (i = 0; i < bw_y; i++) {
        readVectorRescaled(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
      }
      jumpToVector(mc.stego, pos_y);
      for (i = 0; i < bw_y; i++) {
        readVectorRescaled(steg, mc.stego, mc.stego_vectors_down_g + i*dim);
      }
      initDArray(mc.results_c_vs_c_g, SQUARE(mc.cache), tpb, 0.);
      initDArray(mc.results_c_vs_s_g, SQUARE(mc.cache), tpb, 0.);
      initDArray(mc.results_s_vs_s_g, SQUARE(mc.cache), tpb, 0.);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.clean_vectors_right_g, mc.results_c_vs_c_g);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.clean_vectors_down_g, mc.stego_vectors_right_g, mc.results_c_vs_s_g);
      launchGammaKernel(mc, dim, bw_x, bw_y, mc.stego_vectors_down_g, mc.stego_vectors_right_g, mc.results_s_vs_s_g);
      mmdKernel<<<BLOCKS( mc.cache*mc.cache, tpb), tpb>>>(-1.*mc.gamma, mc.results_c_vs_c_g, mc.results_c_vs_s_g, mc.results_s_vs_s_g);
      cudaThreadSynchronize();
      CUBLAS_CALL( cublasGetVector(SQUARE(mc.cache), sizeof(double), mc.results_c_vs_c_g, 1, mc.results, 1));
      for (i = 0ull; i < bw_x; i++) {
        for (j = 0ull; j < bw_y; j++) {
          if (pos_x + i == pos_y + j) continue;
              mmd += mc.results[j + i*mc.cache];
        }
      }
    }
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
  mmd /= (double) (M * (M-1));
  mmd = sqrt(mmd);
  printf("have some mmd: %g [%is]\n", mmd, time(NULL)-start);
  mc.mmd = mmd;
}