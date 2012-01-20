#ifndef STEGOFRAME
#define STEGOFRAME

#include <Stegosaurus.h>
#include "FeatureDock.h"
#include "ConfigDock.h"
#include "GraphDock.h"

#include <QtGui>

const int ranges[] = {152, 48, 34}; // Luma, Chroma DC, Chroma AC

class HistogramWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  HistogramWidget(QWidget *parent = 0, StegoModel *model = 0);
  
  void updateView();
protected:
  int barHeight;
  StegoModel *model;
  void paintEvent(QPaintEvent *event);
  void paintHistogram(int x, int y, double* values, int dim, QPainter *painter);
};


class PairWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  PairWidget(QWidget *parent = 0);
  
  void updateView();
};


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
  QMenu *showMenu;
  QProgressBar *progress;
  QLabel *statusLabel;
  
  HistogramWidget *hw;
  PairWidget *pw;
  
  QAction *openAction;
  QAction *showTableAction;
public slots:
  void openCollection();
};

#endif /* STEGOFRAME */