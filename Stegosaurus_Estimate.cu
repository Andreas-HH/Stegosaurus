#include<Stegosaurus_Estimate.h>
#include<cublas.h>

int init_gpu() {
//   cublasHandle_t handle;
  
  cublasInit();
//   cublasCreate(&handle);
  return 0;
}