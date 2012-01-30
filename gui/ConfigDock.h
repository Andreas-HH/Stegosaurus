#ifndef STEGOCONF
#define STEGOCONF

#include <Stegosaurus.h>
#include "StegoWidgets.h"
#include <QtGui>


class ConfigWidget : public QWidget{
  Q_OBJECT
public:
  ConfigWidget(HistogramWidget *hw = 0, PairWidget *pw = 0, QWidget *parent = 0);
protected:
  HistogramWidget *hw;
  PairWidget *pw;
  QGridLayout *layout;
  QSlider *hdSlider;
  QSlider *scSlider;
  QSlider *bhSlider;
signals:
  void newHorizD(int hd);
  void newScale(int scale);
  void newBarHeight(int bh);
};

class ConfigDock : public QDockWidget {
  Q_OBJECT
public:
  ConfigDock(HistogramWidget *hw = 0, PairWidget *pw = 0, QWidget *parent = 0);
protected:
  QBoxLayout *layout;
  ConfigWidget *panel;
// signals:
//   void newHorizD(int hd);
};

#endif  /*STEGOCONF*/