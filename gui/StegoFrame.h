#ifndef STEGOFRAME
#define STEGOFRAME

#include <Stegosaurus.h>
// #include <StegoClassifier.h>
#include "FeatureDock.h"
#include "ConfigDock.h"
#include "GraphDock.h"
#include "TableDock.h"
#include "StegoWidgets.h"

#include <QtGui>
#include <QtXml>

// const char *slice_types[] = {"P", "B"};

class LoadDialog : public QDialog {
public:
  LoadDialog(QWidget *parent = 0, Qt::WindowFlags f = 0);
  
  QString getMinMethod();
  QString getMaxMethod();
  QString getQPOffset();
  QString getQPRange();
  QString getType();
protected:
  QComboBox *minMethod;
  QComboBox *maxMethod;
  QComboBox *qpOffset;
  QComboBox *qpRange;
  QComboBox *type;
  
  QPushButton *ok;
  QPushButton *cancel;
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
  QMenu *calcMenu;
  QMenu *showMenu;
  QMenu *docMenu;
  QProgressBar *progress;
  QLabel *statusLabel;
  
  HistogramWidget *hw;
  PairWidget *pw;
  
  QFileDialog *fdial;
  QAction *openFeaturesAction;
  QAction *loadFeaturesAction;
  QAction *mmdAction;
  QAction *svmAction;
  QAction *muAction;
  QAction *showTableAction;
  QAction *openDocAction;
  QAction *saveDocAction;
  
  QDomDocument *document;
  QDomElement sets;
  QDomElement mmds;
  QFile *xmlFile;
  LoadDialog *ldial;
public slots:
  void openCollection();
  void calcMMDs();
  void calcMus();
  void saveXML();
  void openXML();
  void classify();
  void loadFeatures();
};

#endif /* STEGOFRAME */