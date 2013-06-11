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

void StegoClassifier::runSVM(double gamma) {
  int sum = 0;
  double C = 100;          // regularization parameter
  bool bias = true;           // use bias/offset parameter
  
  shark::GaussianRbfKernel<> kernel(gamma); // Gaussian kernel
  shark::KernelExpansion<shark::RealVector> ke(&kernel, bias); // (affine) linear function in kernel-induced feature space
  shark::ZeroOneLoss<unsigned int, shark::RealVector> loss; // 0-1 loss

  // define the machine
  shark::CSvmTrainer<shark::RealVector> trainer(&kernel, C);
  shark::ClassificationDataset data(vec, labels);
  shark::ClassificationDataset training, test;
  
  shark::CVFolds<shark::ClassificationDataset> folds = shark::createCVSameSize(data, 10);
  
  training = folds.training(0);
  test     = folds.validation(0);
  std::cout << training.numberOfElements() << std::endl;
  
  for (int i = 0; i < 3; i++) {
    std::cout << "C = " << C << std::endl;
    std::cout << "Algorithm: " << trainer.name() << "\ntraining ..." << std::flush; // Shark algorithms know their names
    trainer.train(ke, training);
    std::cout << "\n  number of iterations: " << trainer.solutionProperties().iterations;
    std::cout << "\n  dual value: " << trainer.solutionProperties().value;
    std::cout << "\n  training time: " << trainer.solutionProperties().seconds << " seconds\ndone." << std::endl;
    
    shark::Data<shark::RealVector> output;  // real-valued output of the machine
    output = ke(training.inputs()); // evaluate on training set
    double train_error = loss.eval(training.labels(), output);
    std::cout << "training error:\t" <<  train_error << std::endl;
    output = ke(test.inputs()); // evaluate on test set
    double test_error = loss.eval(test.labels(), output);
    std::cout << "test error:\t" << test_error << std::endl;
    
    C *= 10.;
    trainer = shark::CSvmTrainer<shark::RealVector>(&kernel, C);
  }
}
