#include <Stegosaurus.h>


int estimateGaussian(stegoContext *steg) {;
  int i;
  int read;
//   double sum;
  int dim = steg->features->dim;
//   double zero = 0.;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
  
//   CUDA_CALL( cudaMemcpy(steg->gpu_c->mu_g, steg->current_feature, dim*sizeof(double), cudaMemcpyHostToDevice));
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), steg->current_feature, 1, steg->gpu_c->vec_g, 1));
//   CUBLAS_CALL( cublasDaxpy(*handle, dim, &zero, steg->gpu_c->vec_g, 1, steg->gpu_c->mu_g, 1));      // init to 0 vector
  for (i = 0; i < steg->features->M; i++) {
    read = readVector(steg->features, steg->current_feature);
    if (read != dim) printf("read something wrong: %i \n", read);
//     sum = 0.;
//     for (j = 0; j < dim; j++)
//       sum += steg->current_feature[j];
//     printf("%f \n", sum);
    CUBLAS_CALL( cublasSetVector(dim, sizeof(double), steg->current_feature, 1, steg->gpu_c->vec_g, 1));
//     CUDA_CALL( cudaMemcpy(steg->gpu_c->vec_g, steg->current_feature, dim*sizeof(double), cudaMemcpyHostToDevice));
    cublasDaxpy(*handle, dim, &(steg->divM), steg->gpu_c->vec_g, 1, steg->gpu_c->mu_g, 1);     // put divM on gpu?
  }
  // print result on screen
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), steg->gpu_c->mu_g, 1, steg->current_feature, 1));
//   CUDA_CALL( cudaMemcpy(steg->current_feature, steg->gpu_c->mu_g, dim*sizeof(double), cudaMemcpyDeviceToHost));
  for (i = 0; i < dim; i++) {
    printf("%f ", steg->current_feature[i]);
//     printf("%f ", steg->divM);
  }
  printf("\n");
  
  return 0;
}
