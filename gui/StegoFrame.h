#ifndef STEGOFRAME
#define STEGOFRAME

#include <Stegosaurus.h>
#include "FeatureDock.h"
#include "ConfigDock.h"
#include "GraphDock.h"
#include "TableDock.h"
#include "StegoWidgets.h"

#include <QtGui>

class StegoFrame : public QMainWindow, public StegoView {
  Q_OBJECT
public:
  StegoFrame(QWidget *parent = 0);
  
  void updateProgress(double p);
protected:
  StegoModel *model;
  FeatureDock *fdock;
  ConfigDock *cdock;
  TableDock *tdock;
  GraphDock *gdock;
  QMenu *fileMenu;
  QMenu *calcMenu;
  QMenu *showMenu;
  QProgressBar *progress;
  QLabel *statusLabel;
  
  HistogramWidget *hw;
  PairWidget *pw;
  
  QFileDialog *fdial;
  QAction *openAction;
  QAction *mmdAction;
  QAction *muAction;
  QAction *showTableAction;
public slots:
  void openCollection();
  void calcMMDs();
  void calcMus();
};

#endif /* STEGOFRAME */