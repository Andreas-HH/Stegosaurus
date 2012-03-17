#include "Stegosaurus.h"


__global__ void calcMMD(int dim, int bw_x, int bw_y, double *down_g, double *right_g, double *results) {
  int i;
  int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
  int idx_y = threadIdx.y + blockIdx.y*blockDim.y;
  int idx_s = threadIdx.y + threadIdx.x*blockDim.y;
  int current_dim = blockDim.x * blockDim.y;
//   int idx_r = idx_x + idx_y*numvec; // one result per block, not per thread!
//   double current_sum = 0.;
  __shared__ double current_sums[1024]; // cuda doesn't allow me to use blockDim.x * blockDim.y here ;_;
  double temp;
  
  if (idx_x < bw_x && idx_y < bw_x) { 
    current_sums[idx_s] = 0.;
    // find || x - y || ^2
    for (i = 0; i < dim; i++) {
      temp = down_g[idx_y*dim + i] - right_g[idx_x*dim + i];
      current_sums[idx_s] += temp * temp;
    }
    // write to global result array
    results[idx_y + bw_y*idx_x] = current_sums[idx_s];
  }
}

__global__ void calcMMD(int dim, int bw_x, int bw_y, double minus_gamma, double *down_g, double *right_g, double *results, bool add) {
  int i;
  int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
  int idx_y = threadIdx.y + blockIdx.y*blockDim.y;
  int idx_s = threadIdx.y + threadIdx.x*blockDim.y;
  int current_dim = blockDim.x * blockDim.y;
//   int idx_r = idx_x + idx_y*numvec; // one result per block, not per thread!
//   double current_sum = 0.;
  __shared__ double current_sums[1024]; // cuda doesn't allow me to use blockDim.x * blockDim.y here ;_;
  double temp;
  bool odd;
  
  if (idx_x < bw_x && idx_y < bw_x) { // they better all pass here, otherwise there will be a deadlock at __syncthreads()
    current_sums[idx_s] = 0.;
    // find || x - y || ^2
    for (i = 0; i < dim; i++) {
      temp = down_g[idx_y*dim + i] - right_g[idx_x*dim + i];
      current_sums[idx_s] += temp * temp;
    }
    // use gamma
    current_sums[idx_s] *= minus_gamma;
    current_sums[idx_s] = exp(current_sums[idx_s]);
    __syncthreads();
    // add the results together
    while (current_dim > 1) {
      odd = current_dim%2 == 1;
      current_dim = current_dim/2;
      if (idx_s < current_dim) {
        current_sums[idx_s] = current_sums[idx_s] + current_sums[current_dim + idx_s];
        if (odd && idx_s == 0) {
	  current_sums[idx_s] = current_sums[idx_s] + current_sums[2*current_dim];
        }
      }
      __syncthreads();
    }
  
    if (idx_s == 0) {
      if (add) results[blockIdx.y + gridDim.y * blockIdx.x] += current_sums[0];
      else     results[blockIdx.y + gridDim.y * blockIdx.x] -= 2*current_sums[0];
    }
  }
//   results[blockIdx.y + gridDim.y * blockIdx.x] = 5;
}

void initMMD(stegoContext *steg, mmdContext& mc) {
  int dim = mc.clean->dim;
  int tpb = steg->gpu_c->threads_per_block;
  
  mc.n = mc.clean->M;
//   mc.cache   = min(mc.n, (int) (steg->gpu_c->doublesOnGPU / (long) dim / 4l)); // 4 input buffers
  mc.kernel_blockwidth = (int) sqrt(steg->gpu_c->threads_per_block);
  mc.cache   = min(mc.n, (int) (sqrt((long) mc.kernel_blockwidth * steg->gpu_c->doublesOnGPU / 3l + SQUARE(2l*(long)mc.kernel_blockwidth*(long)dim / 3l)) - 2l*(long)mc.kernel_blockwidth*(long)dim / 3l)); // 4 input buffers
  mc.kernel_gridwidth = (mc.cache + mc.kernel_blockwidth-1)/mc.kernel_blockwidth;
  printf("cache: %i, kbw: %i \n", mc.cache, mc.kernel_blockwidth);
  
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.clean_vectors_right_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_down_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.stego_vectors_right_g, dim*mc.cache*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_c_vs_c_g, mc.kernel_gridwidth*mc.kernel_gridwidth*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_c_vs_s_g, mc.kernel_gridwidth*mc.kernel_gridwidth*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.results_s_vs_s_g, mc.kernel_gridwidth*mc.kernel_gridwidth*sizeof(double)));
  CUDA_CALL( cudaHostAlloc(&mc.results, mc.kernel_gridwidth*mc.kernel_gridwidth*sizeof(double), cudaHostAllocDefault));
//   CUDA_CALL( cudaHostAlloc(&mc.stego_vectors, dim*mc.cache*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mc.v_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&mc.v2_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&mc.temp_g, dim*sizeof(double)));
  
//   initDArray<<<BLOCKS(mc.kernel_gridwidth*mc.kernel_gridwidth, tpb), tpb>>>(mc.results_c_vs_c_g, mc.kernel_gridwidth*mc.kernel_gridwidth, 1.);
  
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
  CUDA_CALL( cudaFree(mc.results_c_vs_c_g));
  CUDA_CALL( cudaFree(mc.results_c_vs_s_g));
  CUDA_CALL( cudaFree(mc.results_s_vs_s_g));
  CUDA_CALL( cudaFree(mc.v_g));
  CUDA_CALL( cudaFree(mc.temp_g));
  CUDA_CALL( cudaFreeHost(mc.results));
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
  int gridwidth, gridheight;
//   dim3 grid;
//   dim3 block;
  long pos_x, pos_y;
  int tpb = steg->gpu_c->threads_per_block;
//   int cache = mc->cache;
//   int blockwidth_g, pos_g, cache_g = mc->cache_g;
//   int lessThanOne = 0, greaterThanOne = 0;
//   double *vectors;
  int dim = mc.clean->dim;
  int M = mc.n;
  double mmd = 0.;
  double temp = 0.;
  
  for (pos_x = 0l; pos_x < (long) M; pos_x += mc.cache) {
    bw_x = min(mc.cache, M-(int)pos_x);
//     printf("clean [%i] \n", pos_x);
    jumpToVector(mc.clean, pos_x);
//     printf("stego [%i] \n", pos_x);
    jumpToVector(mc.stego, pos_x);
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.clean, mc.clean_vectors_down_g + i*dim);
    } // better to one after the other, might be quicker
    for (i = 0; i < bw_x; i++) {
      readVectorRescaled(steg, mc.stego, mc.stego_vectors_down_g + i*dim);
    }
    for (pos_y = 0l; pos_y < (long) M; pos_y += mc.cache) {
      bw_y = min(mc.cache, M-(int)pos_y);
//       printf(" clean [%i] \n", pos_y);
      jumpToVector(mc.clean, pos_y);
//       printf(" stego [%i] \n", pos_y);
      jumpToVector(mc.stego, pos_y);
      for (i = 0; i < bw_x; i++) {
        readVectorRescaled(steg, mc.clean, mc.clean_vectors_right_g + i*dim);
      }
      for (i = 0; i < bw_x; i++) {
        readVectorRescaled(steg, mc.stego, mc.stego_vectors_right_g + i*dim);
      }
//       for (j = 0; j < bw_x; j++) {
// 	for (k = 0; k < bw_y; k++) {
// 	  if (pos_x+j == pos_y + k)  continue;
// 	  mmd += applyKernel(steg, mc.gamma, dim, mc.clean_vectors_down_g+j*dim, mc.clean_vectors_right_g+k*dim, mc.temp_g);
// 	  mmd += applyKernel(steg, mc.gamma, dim, mc.stego_vectors_down_g+j*dim, mc.stego_vectors_right_g+k*dim, mc.temp_g);
// 	  mmd -= 2.*applyKernel(steg, mc.gamma, dim, mc.clean_vectors_down_g+j*dim, mc.stego_vectors_right_g+k*dim, mc.temp_g);
// 	}
//       }
      gridwidth = (bw_x + mc.kernel_blockwidth-1)/mc.kernel_blockwidth;
      gridheight = (bw_y + mc.kernel_blockwidth-1)/mc.kernel_blockwidth;
      initDArray<<<BLOCKS(gridwidth*gridheight, tpb), tpb>>>(mc.results_c_vs_c_g, mc.kernel_gridwidth*mc.kernel_gridwidth, 0.);
      calcMMD<<<(gridwidth, gridheight), (mc.kernel_blockwidth, mc.kernel_blockwidth)>>>(dim, bw_x, bw_y, -1.*0.2, mc.clean_vectors_down_g, mc.clean_vectors_right_g, mc.results_c_vs_c_g, true);
      calcMMD<<<(gridwidth, gridheight), (mc.kernel_blockwidth, mc.kernel_blockwidth)>>>(dim, bw_x, bw_y, -1.*0.2, mc.stego_vectors_down_g, mc.stego_vectors_right_g, mc.results_c_vs_c_g, true);
      calcMMD<<<(gridwidth, gridheight), (mc.kernel_blockwidth, mc.kernel_blockwidth)>>>(dim, bw_x, bw_y, -1.*0.2, mc.clean_vectors_down_g, mc.stego_vectors_right_g, mc.results_c_vs_c_g, false);
//       CUBLAS_CALL( cublasGetMatrix(gridheight, gridwidth, sizeof(double), mc.results_c_vs_c_g, mc.kernel_gridwidth, mc.results, mc.kernel_gridwidth));
      CUBLAS_CALL( cublasGetVector(SQUARE(mc.kernel_gridwidth), sizeof(double), mc.results_c_vs_c_g, 1, mc.results, 1));
      temp = 0.;
      printf("results[10]: %g \n", mc.results[10]);
      for (j = 0; j < gridwidth; j++) {
	for (k = 0; k < gridheight; k++) {
	  if (pos_x+j == pos_y + k)  continue;
	  temp += mc.results[j*mc.kernel_gridwidth + k];
//           printf("%g ", mc.results[j*gridwidth + k]);
	  mmd += mc.results[j*gridwidth + k];
	}
      }
      printf("delta mmd: %f \n", temp);
    }
//     stegoRewind(mc.stego);
  }
  stegoRewind(mc.clean);
  stegoRewind(mc.stego);
  printf("have some mmd: %f \n", mmd);
  mc.mmd = 1./((double) M*(M-1)) * mmd;
}