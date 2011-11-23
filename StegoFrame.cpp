#include<Stegosaurus_Estimate.h>
#include<Stegosaurus_KL.h>

#include<QApplication>
#include<QLabel>


int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  QLabel *label = new QLabel("Stegosaurus");
  label->show();
  estimateGaussian("data/p_histograms.fv");
  return app.exec();
}