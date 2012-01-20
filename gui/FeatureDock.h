#ifndef STEGOFDOCK
#define STEGOFDOCK

#include <Stegosaurus.h>
#include <QtGui>


class FeatureWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  FeatureWidget(QWidget* parent = 0, StegoModel *model = 0);
  
//   void updateView();
  void updateCollection();
//   void updateProgress(double p);
//   void setModel(StegoModel *model);
protected:
  QBoxLayout *layout;
  StegoModel *model;
//   QProgressBar *progress;
  FeatureCollection *collection;
};

class FeatureDock : public QDockWidget {
  Q_OBJECT
public:
  FeatureDock(QWidget *parent = 0, StegoModel *model = 0);
};

#endif  /*STEGOFDOCK*/