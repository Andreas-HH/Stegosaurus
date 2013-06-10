#include "StegoClassifier.h"

StegoClassifier::StegoClassifier(int dim) {
  this->dim = dim;
}

void StegoClassifier::addVector(double *data) {
  shark::RealVector rv(dim);
  
  for (int i = 0; i < dim; i++) {
    rv(i) = data[i];
//     std::cout << rv(i) << " ";
  }
//   std::cout << std::endl;
  
  vec.push_back(rv);
}

void StegoClassifier::addCleanVector(double *data) {
  addVector(data);
  labels.push_back(0);
}

void StegoClassifier::addStegoVector(double *data) {
  addVector(data);
  labels.push_back(1);
}
