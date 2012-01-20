#ifndef STEGOSAURUS
#define STEGOSAURUS

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cublas_v2.h>
#include <cuda.h>

#define CUBLAS_CALL(x) switch (x) {       \
    case CUBLAS_STATUS_SUCCESS:           \
      break;                              \
    case CUBLAS_STATUS_NOT_INITIALIZED:   \
      printf("not init \n");              \
      break;                              \
    case CUBLAS_STATUS_INVALID_VALUE:     \
      printf("invalid value \n");         \
      break;                              \
    case CUBLAS_STATUS_MAPPING_ERROR:     \
      printf("Mapping error \n");         \
      break;                              \
    default:                              \
      printf("Something else \n");        \
      break;                              \
  }
  
#define CUDA_CALL(x) switch (x) {                        \
    case cudaSuccess:                                    \
      break;                                             \
    case cudaErrorInvalidValue:                          \
      printf("cuda: invalid value! \n");                 \
      break;                                             \
    case cudaErrorInvalidDevicePointer:                  \
      printf("cuda: invalid dev pointer! \n");           \
      break;                                             \
    case cudaErrorInvalidMemcpyDirection:                \
      printf("cuda: invalid memcopy direction! \n");     \
      break;                                             \
    default:                                             \
      printf("cuda: Something else. \n");                \
      break;                                             \
  }

#define BLOCKS(x,tbp) ((x)+(tbp)-1)/(tpb)

const int QP_RANGE = 8;


typedef struct myGaussian {
  int dim;
  double *mu;
  double *sigma;
  double *sigma_inverse;
  double *qr_diag;
  double *qr;
} myGaussian;

typedef struct featureSet {
  int dim;
  int M;
  int gpu_matrix_width;                      // max number of matrix column in gpu memory
  FILE *file;
  char *name;
//   double *current_feature;                   // created using cudaHostAlloc (important!) 
  double *max_vec;                           // vector of maximum entries
  double *qp_vec;
  
//   double *qp_g;
//   double *mu_g;
//   double *sigma_g;
//   double *vec_g;
//   double *max_g;                             // maybe a waste of gpu memory when it comes to computing KL
  
  double kl_div;
  double divM;    
  myGaussian* gauss;
} featureSet;

typedef struct gpuContext {
  int threads_per_block;
  cublasHandle_t handle;
} gpuContext;

typedef struct stegoContext {
  gpuContext *gpu_c;                         // might contain multiple cublas handles later (or multiple gpucontexts, we will see...)
  featureSet *features;                         // just have an array of handles here?
//   myGaussian *clean_gauss;                   // shared among stegoContexts, only read though! if multiple contexts are loaded at all
//   myGaussian *stego_gauss;                   // doesn't store inverse or QR.
} stegoContext;


class StegoView {
public:
  virtual void updateView() { };
  virtual void updateProgress(double p) { };
  virtual void updateCollection() { };
};

class FeatureCollection {
public:
  FeatureCollection(char *path);
  ~FeatureCollection();
  
  void rewind();
  void setSelected();
  int isSelected();
  int getNumSets();
  int getCurrentSet();
  featureSet* nextFeatureSet();
  int hasNext();
  featureSet **collection;                   // Array of featureSet's
protected:
  int num_sets;
  int current_set;
  int selected;
};

class StegoModel {
public:
  StegoModel();
  ~StegoModel();
  
  void addView(StegoView *view);
  void openCollection(char *path);
  void estimateMus();
  
  int getDimension();
  double* getMaxVector();                    // Vector of maximum elements
  double* getMuVector();                     // can be a particular feature vector or mu
  double* getQPHist();
  double* getSigma();
  FeatureCollection* getCollection();
protected:
  int current_view;
//   double progress;
  StegoView *views[25];                      // maybe a linked list or some other container is better here
  FeatureCollection *fcol;                     // or some signal/slot system
  stegoContext *steg;
  void modelChanged();                       // asks all views to update themselves
  void progressChanged(double p);
  void collectionChanged();
};


// void storeGaussian(char *path, myGaussian *gauss);
featureSet* openFeatureSet(char *path);
int closeFeatureSet(featureSet *set);
int readVector(featureSet *set, double *vec);
void stegoRewind(featureSet *set);
void pathConcat(char* a, char* b, char *result);

extern "C" { 
  stegoContext* init_stego();
  gpuContext* init_gpu();
  void close_stego(stegoContext *steg);
  void close_gpu(gpuContext *gp);
  
  int estimateMu(stegoContext *steg);
  int estimateSigma(stegoContext *steg);
//   int computeQPHistogram(stegoContext *steg, double *mu_g, int qp_range, double *result);
}

#endif /* STEGOSAURUS */
