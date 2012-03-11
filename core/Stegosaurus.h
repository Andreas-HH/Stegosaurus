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
#include <queue>
// #include <unordered_map>

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

#define BLOCKS(x,tpb) ((x)+(tpb)-1)/(tpb)

using namespace std;

const int QP_RANGE = 20;
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
  int video_bitrate;
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
  int id;                                    // helps to find a set
  int gpu_matrix_width;                      // max number of matrix column in gpu memory
  int num_files;
  int current_file;
  long *vsPerFile;
  FILE **files;
  char *name;
  long dataOffset;

  double *vec;
  double *vec_g;
  double *ones_g;
  double *max_vec;                           // vector of maximum entries
  double *qp_vec;
  double rate;
  double divM;    
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
  int n;
  int cache;
  double gamma;
  double mmd;
  featureSet *clean;
  featureSet *stego;
  double *clean_vectors_down_g;
  double *clean_vectors_right_g;
  double *stego_vectors_down_g;
  double *stego_vectors_right_g;
//   double *clean_vectors;
  double *v_g;
//   double *v2_g;
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
  gpuContext *gpu_c;                         // might contain multiple cublas handles later (or multiple gpucontexts, we will see...)
  featureSet *features;                      // just have an array of handles here?
  long doublesInRAM;
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
  
//   class Iterator {
//   public:
//     Iterator(StegoModel *model);
// //     void increaseLevel();
// //     void decreaseLevel();
//     int getLevel0();
//     int getLevel1();
//     int getLevel2();
//     int getLevel3();
// //     FeatureCollection::Iterator *nextCollectionIterator();
//     bool hasNext();
//   protected:
//     int level;
//     int level1, level3;
//     map< int, vector< map< int, vector< FeatureCollection* > > > >::iterator *level0iter;
//     vector< map< int, vector< FeatureCollection* > > >::iterator             *level1iter;
//     map< int, vector< FeatureCollection* > >::iterator                       *level2iter;
//     vector< FeatureCollection* >::iterator                                   *level3iter;
//     void next();
//   };
  
  void addView(StegoView *view);
  void openDirectory(const char* path);
  void estimateMus();
  double doMMD(featureSet* clean, featureSet* stego);
  void setFeatures(featureSet *set);
  featureSet* getCleanSet(); // there must be some class-frontend eventually to allow DB!!!
  
  int getDimension();
  double* getMaxVector();                    // Vector of maximum elements
  double* getMuVector();                     // can be a particular feature vector or mu
  double* getQPHist();
  double* getSigma();
  double *getDiag();
  int** getRanges();
//   Iterator* getIterator();
  FeatureCollection::Iterator* getFeatureIterator(int video_birate, int pair, int method, int accept);
  FeatureCollection* getCollection();
protected:
//   int current_view;
  int **ranges;                              //  [block][row], row depends on if we have pair or hist features,, the first added feature file will determine this!
  list< StegoView* > views;
  FeatureCollection *collections[10][8];           // [method][accept] ... [0] = clean, [1] = plusminus1, ... 
//   map< int, vector< map< int, vector< FeatureCollection* > > > > collections; // video_bitrate -> hist/pair -> method -> accept
  stegoContext *steg;
  featureSet *cleanSet;
  mmdContext *mc;
  void modelChanged();                       // asks all views to update themselves
  void progressChanged(double p);
  void collectionChanged();
};


// void storeGaussian(char *path, myGaussian *gauss);
int readHeader(FILE *file, featureHeader *header);
int closeFeatureSet(featureSet *set);
int jumpToVector(featureSet *set, long vecnum);
int newFeatureFile(featureSet *set, const char *path);
void stegoRewind(featureSet *set);
void pathConcat(const char* a, const char* b, char *result);

extern "C" { 
  stegoContext* init_stego();
  gpuContext* init_gpu();
  void close_stego(stegoContext *steg);
  void close_gpu(gpuContext *gp);
  int estimateMu(stegoContext *steg);
  void printPointerInfo(void *ptr);
  featureSet* openFeatureSet(const char* path);
  int readVector(stegoContext* steg, featureSet* set, double* vec);
  int readVectorL2(stegoContext* steg, featureSet* set, double* vec);
  
  void initMMD(stegoContext* steg, mmdContext& mc);
  void closeMMD(mmdContext& mc);
  void estimateGamma(stegoContext* steg, mmdContext& mc);
  void estimateMMD(stegoContext* steg, mmdContext& mc);
  double applyKernel(stegoContext* steg, double gamma, int dim, double* v1_g, double* v2_g, double* temp_g);
  
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
