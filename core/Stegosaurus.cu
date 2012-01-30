#include <Stegosaurus.h>


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

//   CUDA_CALL( cudaFreeHost(set->current_feature));
//   CUDA_CALL( cudaFreeHost(set->max_vec));
//   CUDA_CALL( cudaFreeHost(set->qp_vec));
//   CUDA_CALL( cudaFree(set->mu_g));
//   CUDA_CALL( cudaFree(set->qp_g));
//   CUDA_CALL( cudaFree(set->vec_g));
//   CUDA_CALL( cudaFree(set->max_g));
//   if (set->gauss->mu != NULL)             free (set->gauss->mu);
  if (set->gauss->mu != NULL)             CUDA_CALL( cudaFree(set->gauss->mu));
  if (set->max_vec != NULL)               CUDA_CALL( cudaFree(set->max_vec));
  if (set->qp_vec != NULL)                CUDA_CALL( cudaFree(set->qp_vec));
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
//   gp->qr_cache = 10;
//   printf("%i \n", gp->threads_per_block);

  return gp;
}

stegoContext* init_stego() {
  stegoContext *steg = (stegoContext*) malloc(sizeof(stegoContext));
  steg->gpu_c = init_gpu();
  steg->features = NULL;

  return steg;
}


void close_stego(stegoContext *steg) {
  close_gpu(steg->gpu_c);
  closeFeatureSet(steg->features);
}

void close_gpu(gpuContext* gp) {
  cublasDestroy(gp->handle);
}
