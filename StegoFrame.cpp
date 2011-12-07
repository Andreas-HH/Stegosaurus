#include <Stegosaurus.h>

#include <QApplication>
#include <QLabel>


int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  QLabel *label = new QLabel("Stegosaurus");
  label->show();
  
  stegoContext *steg =  init_stego("p_histograms.fv");

  estimateGaussian(steg);
  
  close_stego(steg);
  
  return app.exec();
}
