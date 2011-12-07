#ifndef STEGOSAURUS
#define STEGOSAURUS

#include <stdio.h>
#include <stdlib.h>
#include <cublas_v2.h>
#include <cuda.h>

#define CUDA_CALL(x) if (x != cudaSuccess) printf("not a cudaSuccess! \n")
#define CUBLAS_CALL(x) if (x != CUBLAS_STATUS_SUCCESS) printf("not a cublasSuccess! \n")



typedef struct featureSet {
  int dim;
  int M;
  FILE *file;
  double *current_feature;
} featureSet;


typedef struct myGaussian {
  int dim;
  double *mu;
  double *sigma;
  double *sigma_inverse;
  double *qr_diag;
  double *qr;
} myGaussian;


typedef struct gpuContext {
  cublasHandle_t handle;
  double *mu_g;
  double *sigma_g;
  double *vec_g;    // feature vector on gpu
//   cublasPointerMode_t pm;
} gpuContext;


// one stegocontext per stego-featureset
typedef struct stegoContext {
  gpuContext *gpu_c;       // might contain multiple cublas handles later (or multiple gpucontexts, we will see...)
  featureSet *features;
  myGaussian *clean_gauss; // shared among stegoContexts, only read though! if multiple contexts are loaded at all
  myGaussian *stego_gauss; // doesn't store inverse or QR.
  
  double kl_div;
  double divM;  // 1/M where M is from the featureset
  double *current_feature; // created using cudaHostAlloc (important!) 
} stegoContext;



void            storeGaussian(char *path, myGaussian *gauss);

featureSet*     openFeatureSet(char *path);
int             closeFeatureSet(featureSet *set);
int             readVector(featureSet *set, double *vec);
void            stegoRewind(featureSet *set);


extern "C" { 
  stegoContext* init_stego(char *features_path);
  gpuContext*   init_gpu(int dim);
  
  void          close_stego(stegoContext *steg);
  void          close_gpu(gpuContext *gp);
  
  int           estimateGaussian(stegoContext *steg);
}

#endif /* STEGOSAURUS */