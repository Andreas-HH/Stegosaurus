#ifndef STEGO_SVM
#define STEGO_SVM

#include "vector"
#include <shark/Data/Dataset.h>
#include <shark/Algorithms/Trainers/CSvmTrainer.h>
#include <shark/Models/Kernels/GaussianRbfKernel.h>
#include <shark/ObjectiveFunctions/Loss/ZeroOneLoss.h>
#include <shark/Data/CVDatasetTools.h>

class StegoClassifier {
public:
  StegoClassifier(int dim);
  
  void addCleanVector(double *data);
  void addStegoVector(double *data);
  
  void runSVM(double gamma);
protected:
  int dim;
  
  std::vector<shark::RealVector>  vec;
  std::vector<unsigned int>       labels;
  
  void addVector(double *data);
};


#endif /* STEGO_SVM */
