#ifndef STEGOSAURUS
#define STEGOSAURUS

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cublas_v2.h>
#include <cuda.h>
#include <vector>
#include <list>
#include <map>

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

using namespace std;

const int QP_RANGE = 8;
const int MAX_FILES = 20;
const int num_coefs[] = {16, 4, 15};

const double BIN_WIDTH = 0.001;

typedef struct myGaussian {
  int dim;
  double *mu;
  double *sigma;
  double *sigma_inverse;
  double *qr_diag;
  double *qr;
} myGaussian;

typedef struct featureHeader {
  char pair;
  char slice_type;
  char method;
  char using_rate;
  double rate;
  char accept;
  char qp_offset;
  char qp_range;
  unsigned char ranges[3][16];
} featureHeader;

typedef struct featureSet {
  featureHeader *header;
  int dim;
  int M;
  int gpu_matrix_width;                      // max number of matrix column in gpu memory
  int num_files;
  int current_file;
  FILE **files;
  char *name;

  double *max_vec;                           // vector of maximum entries
  double *qp_vec;
  double kl_div;
  double divM;    
  myGaussian* gauss;
} featureSet;

typedef struct gpuContext {
  int threads_per_block;
  int num_streams;
  cublasHandle_t handle;
} gpuContext;

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
  
} mmdContext;

typedef struct stegoContext {
  gpuContext *gpu_c;                         // might contain multiple cublas handles later (or multiple gpucontexts, we will see...)
  featureSet *features;                      // just have an array of handles here?
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
  int getCurrentSet();
  int addFeatureFile(const char *path, featureHeader *header);
  featureSet* getFeatureSet(int index);
  featureSet* nextFeatureSet();
  FeatureCollection::Iterator* iterator();
  int hasNext();
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
  
  void addView(StegoView *view);
  void openDirectory(const char* path);
  void estimateMus();
  
  int getDimension();
  double* getMaxVector();                    // Vector of maximum elements
  double* getMuVector();                     // can be a particular feature vector or mu
  double* getQPHist();
  double* getSigma();
  double *getDiag();
  int** getRanges();
  FeatureCollection* getCollection();
protected:
//   int current_view;
  int **ranges;                              //  [block][row], row depends on if we have pair or hist features,, the first added feature file will determine this!
  list< StegoView* > views;
  FeatureCollection *collections[10][8];           // [method][accept] ... [0] = clean, [1] = plusminus1, ... 
  stegoContext *steg;
  void modelChanged();                       // asks all views to update themselves
  void progressChanged(double p);
  void collectionChanged();
};


// void storeGaussian(char *path, myGaussian *gauss);
int readHeader(FILE *file, featureHeader *header);
featureSet* openFeatureSet(const char* path);
int closeFeatureSet(featureSet *set);
int newFeatureFile(featureSet *set, const char *path);
int readVector(featureSet *set, double *vec);
void stegoRewind(featureSet *set);
void pathConcat(const char* a, const char* b, char *result);

extern "C" { 
  stegoContext* init_stego();
  gpuContext* init_gpu();
  void close_stego(stegoContext *steg);
  void close_gpu(gpuContext *gp);
  int estimateMu(stegoContext *steg);
  void printPointerInfo(void *ptr);
  
  int estimateSigma(stegoContext *steg);
  void constructQ(double *q, double *diag, double norm, int count);
  klContext* initKLContext(stegoContext* steg);
  void closeKLContext(klContext *klc);
  void applyQ(cublasHandle_t handle, int dim, double *q, double *vec, double *tmp, double *dotp, double *mintwo);
  void applyQs(klContext* klc, double* qs_g, double* vs_g, int numqs, int numvs, int cdim);
  void qrHouseholder(stegoContext *steg);
  void qrUnitTest();
}

#endif /* STEGOSAURUS */
