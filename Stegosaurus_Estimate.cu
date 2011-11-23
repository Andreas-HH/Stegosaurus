#include<Stegosaurus_Estimate.h>
#include<cublas.h>

int init_gpu() {
//   cublasHandle_t handle;
  
  cublasInit();
  printf("could do something useful on the gpu now...\n");
  cublasShutdown();
//   cublasCreate(&handle);
  return 0;
}
