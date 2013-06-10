#ifndef STEGOSAURUS
#define STEGOSAURUS

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
// #include <strstream>
#include <cublas_v2.h>
#include <cuda.h>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <utility>
#include <string>
#include <stdint.h>
// #include <unordered_map>

#define SQUARE(x) (x)*(x)

#define CUBLAS_CALL(x) switch (x) {       \
    case CUBLAS_STATUS_SUCCESS:           \
      break;                              \
    case CUBLAS_STATUS_NOT_INITIALIZED:   \
      printf("cublas: not init \n");              \
      break;                              \
    case CUBLAS_STATUS_INVALID_VALUE:     \
      printf("cublas: invalid value \n");         \
      break;                              \
    case CUBLAS_STATUS_MAPPING_ERROR:     \
      printf("cublas: Mapping error \n");         \
      break;                              \
    default:                              \
      printf("cublas: Something else \n");        \
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

#define BLOCKS(x,tpb) ((x)+(tpb)-1)/(tpb)

using namespace std;

const int QP_RANGE = 20;
const int MAX_FILES = 20;
const int num_coefs[] = {16, 4, 15};

// const double BIN_WIDTH = 0.001;
typedef uint32_t store_elem;

typedef struct myGaussian {
  int dim;
  double *mu;
  double *sigma;
  double *sigma_inverse;
  double *qr_diag;
  double *qr;
} myGaussian;

typedef struct featureHeader {
  int video_bitrate;
//   char pair;   // not necessary
  char slice_type;
  char method;
//   char using_rate;
  double prob;
  char accept;
  char qp_offset;
  char qp_range;
  unsigned char ranges[3][16];
} featureHeader;

typedef struct featureSet {
  featureHeader *header;
  int dim;
  int dim_file;
  int hist_dim; // each for single qp
  int pair_dim;
  int uvsv_dim;
  int masked_dim;
  uint64_t M;
  int id;                                    // helps to find a set
  int gpu_matrix_width;                      // max number of matrix column in gpu memory
  int num_files;
  int current_file;
  uint64_t *vsPerFile;
  FILE **files;
  char **paths;
  char *name;
  long dataOffset;

  store_elem *counts;                               // want to make this long as soon as compression works
  double *vec;
  double *vec_g;
  double *ones_g;
  double *max_g;
  double *min_g;                             // usually is identically zero, but better be sure
  double *mu_g;
  double *mu_vec;
  double *var_g;
  uint64_t *mask_counts;
  int *mask_vec;
  double *max_vec;                           // vector of maximum entries
  double *min_vec;
  double *qp_vec;
  double prob;
  double divM;
  double mmd;
  myGaussian* gauss;
} featureSet;

typedef struct klContext {
  int dim;
  int nstr;

  cudaStream_t *streams;
  double **tmp_g;
  double **vec_g;
  double **dotp_g;
  cublasHandle_t *handle;
} klContext;

typedef struct mmdContext {
  uint64_t n;
  uint64_t cache;
  int kernel_blockwidth;
  int kernel_gridwidth;
  double gamma;
  double mmd;
  featureSet *clean;
  featureSet *stego;
  double *clean_vectors_down_g;
  double *clean_vectors_right_g;
  double *stego_vectors_down_g;
  double *stego_vectors_right_g;
  double *results_c_vs_c_g;
  double *results_c_vs_s_g;
  double *results_s_vs_s_g;
  double *results;
  double *v_g;
  double *temp_g;
  double *vectors_g;
} mmdContext;

typedef struct gpuContext {
  int threads_per_block;
  int num_streams;
  long doublesOnGPU;
  cublasHandle_t handle;
} gpuContext;

typedef struct stegoContext {
  gpuContext *gpu_c;                    // might contain multiple cublas handles later (or multiple gpucontexts, we will see...)
  featureSet *features;                 // just have an array of handles here?
  long doublesInRAM;
} stegoContext;


class StegoView {
public:
  virtual void updateView() { };
  virtual void updateProgress(double p) { };
  virtual void updateCollection() { };
};

class FeatureCollection {
public:
  FeatureCollection(featureHeader *h);
  ~FeatureCollection();
  
  class Iterator {
  public:
    Iterator(FeatureCollection *f);
    bool hasNext();
    featureSet* next();
  protected:
    FeatureCollection *fc;
    map< int, featureSet* >::iterator iter;
  };
  
  int getNumSets();
  int addFeatureFile(const char *path, featureHeader *header, stegoContext *steg, featureSet *cleanSet);
  featureSet* getFeatureSet(int index);
  FeatureCollection::Iterator* iterator();
protected:
  map < int, featureSet* > collection;
  featureHeader header;
  int num_sets;
  int current_set;
  int selected;
};

class StegoModel {
public:
  StegoModel();
  ~StegoModel();
  
  class Iterator {
  public:
    Iterator(StegoModel* m);
    bool hasNext();
    pair< int, int > getX();
    FeatureCollection* next();
  protected:
    StegoModel *model;
    pair< int, int > x;
    map< pair< int, int >, FeatureCollection* >::iterator iter;
  };
  
  void runClassifier();
 
  void addView(StegoView *view);
  void openDirectory(const char *path);
  int openFile(const char* path, int i, int num_sets, featureHeader& header);
  void estimateMus();
  void doMMDs();
  double doMMD(featureSet* clean, featureSet* stego);
  void setFeatures(featureSet *set);
  featureSet* getCleanSet(); // there must be some class-frontend eventually to allow DB!!!
  
  int getDimension();
  int getHistDim();
  int getPairDim();
  int getUvsVDim();
  int getSigmaDim();
  int getQPRange();
  double* getMaxVector();                    // Vector of maximum elements
  double* getMuVector();                     // can be a particular feature vector or mu
  double* getQPHist();
  double* getSigma();
  double *getDiag();
  int** getRanges();
  StegoModel::Iterator* iterator();
  
  void collectionChanged(); // rebuilds the collection tree in the GUI
protected:
  int **ranges;                             
  list< StegoView* > views;
  map< pair< int, int >, FeatureCollection* > collections;
  set<string> *seenPaths;
  stegoContext *steg;
  featureSet *cleanSet;
  mmdContext *mc;
  void modelChanged();                      // asks all views to update themselves
  void progressChanged(double p);           // aks all views to update their progress
};


// void storeGaussian(char *path, myGaussian *gauss);
int readHeader(FILE *file, featureHeader *header);
void startAction(featureSet *set);  // opens the files of the given featureSet
void endAction(featureSet *set);    // closes the corresponding files
int closeFeatureSet(featureSet *set); // frees all resources related to the given featreuSet
int jumpToVector(featureSet* set, uint64_t vecnum);
void stegoRewind(featureSet *set);
void pathConcat(const char* a, const char* b, char *result);

// normalization kernels
__global__ void initDArrayKernel(double *m, int dim, double val);
__global__ void finishMax(int dim, double* min, double* max);
__global__ void compareMax(int dim, double *current_max, double *new_features);
__global__ void compareMin(int dim, double *current_min, double *new_features);
__global__ void rescaleKernel(int dim, double *vec_g, double *min_g, double *max_g);
__global__ void varianceKernel(double divM, double *vec_g, double *mu_g, double *var_g, int dim);
__global__ void normalizeKernel(double* vec_g, double* mu_g, double* var_g, int dim);

// MMD kernels
__global__ void gammaKernel(int dim, uint64_t cache, int offset, int steps, uint64_t bw_x, uint64_t bw_y, double* down_g, double* right_g, double* results);
__global__ void mmdKernel(double minus_gamma, double *cvc_g, double *cvs_g, double *svs_g);

// QR factorization kernels
__global__ void subtract(int dim, double *vec, double *mu);
__global__ void constructQ(double *q, double *diag, double norm, int count);

extern "C" { 
  void initDArray(double* m, int dim, int tpb, double val);
  stegoContext* init_stego();
  gpuContext* init_gpu();
  void close_stego(stegoContext *steg);
  void close_gpu(gpuContext *gp);
  int estimateMu(stegoContext *steg);
  void printPointerInfo(void *ptr);
  featureSet* openFeatureSet(const char* path, stegoContext *steg);
  int newFeatureFile(stegoContext* steg, featureSet* set, const char* path);
  int estimateScalingParameters(stegoContext* steg, featureSet* set);
//   int readCountVector(double *data, int *cache, int dim, FILE *file);
  int readCounts(featureSet* set);
  int readVectorRescaled(stegoContext *steg, featureSet *set, double *vec_g);
  int readVectorNormalized(stegoContext *steg, featureSet *set, double *vec_g);
  int readVectorL1D(stegoContext *steg, featureSet *set, double *vec_g);
  int readVectorL2(stegoContext *steg, featureSet *set, double *vec_g);
  void scaleL1D(stegoContext* steg, int dim, double* vec, double* vec_g, double* ones_g);
  
  void initMMD(stegoContext* steg, mmdContext& mc);
  void closeMMD(mmdContext& mc);
  void estimateGamma(stegoContext* steg, mmdContext& mc);
  void launchGammaKernel(mmdContext& mc, int dim, uint64_t bw_x, uint64_t bw_y, double* down_g, double* right_g, double* results_g);
  void estimateMMD(stegoContext* steg, mmdContext& mc);
//   double applyKernel(stegoContext* steg, double gamma, int dim, double* v1_g, double* v2_g, double* temp_g);
  
  int estimateSigma(stegoContext *steg);
  klContext* initKLContext(stegoContext* steg);
  void closeKLContext(klContext *klc);
  void applyQ(cublasHandle_t handle, int dim, double *q, double *vec, double *tmp, double *dotp, double *mintwo);
  void applyQs(klContext* klc, double* qs_g, double* vs_g, int numqs, int numvs, int cdim);
  void qrHouseholder(stegoContext *steg);
  void qrUnitTest();
}

#endif /* STEGOSAURUS */
