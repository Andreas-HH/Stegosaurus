#ifndef STEGO_SVM
#define STEGO_SVM

#include "vector"
#include "shark/Data/Dataset.h"

class StegoClassifier {
public:
  StegoClassifier(int dim);
  
  void addCleanVector(double *data);
  void addStegoVector(double *data);
  
  void validate(double gamma);
protected:
  int dim;
  
  std::vector<shark::RealVector> vec;
  std::vector<unsigned int>      labels;
  shark::ClassificationDataset data;
  
  void addVector(double *data);
};


#endif /* STEGOSAURUS */
