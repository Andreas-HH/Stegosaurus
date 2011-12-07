#include <Stegosaurus.h>


gpuContext* init_gpu(int dim) {
  gpuContext *gp = (gpuContext*) malloc(sizeof(gpuContext));
//   cublasHandle_t handle;
//   cudaError_t cerr;

  cublasCreate(&(gp->handle));
  CUDA_CALL( cudaMalloc(&(gp->mu_g), dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&(gp->vec_g), dim*sizeof(double)));

//   if (cerr != cudaSuccess) printf("not a cudaSuccess! \n");
  
//   printf("could do something useful on the gpu now...\n");
  return gp;
}

stegoContext* init_stego(char *features_path) {
  int i;
  stegoContext *steg = (stegoContext*) malloc(sizeof(stegoContext));
  
  steg->features = openFeatureSet(features_path);
  steg->gpu_c = init_gpu(steg->features->dim);
  CUDA_CALL( cudaHostAlloc(&(steg->current_feature), steg->features->dim*sizeof(double), cudaHostAllocDefault));
  
  for (i = 0; i < steg->features->dim; i++)
    steg->current_feature[i] = 0.;
  
  steg->divM = 1./((double) steg->features->M);
  return steg;
}


void close_stego(stegoContext *steg) {
  close_gpu(steg->gpu_c);
  closeFeatureSet(steg->features);
  CUDA_CALL( cudaFreeHost(steg->current_feature));
}

void close_gpu(gpuContext* gp) {
  cublasDestroy(gp->handle);
  CUDA_CALL( cudaFree(gp->mu_g));
  CUDA_CALL( cudaFree(gp->vec_g));
}
