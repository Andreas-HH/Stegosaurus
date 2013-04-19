#ifndef STEGOWIDGETS
#define STEGOWIDGETS

#include <Stegosaurus.h>
#include <QtGui>
#include <cmath>


// const int hist_blocktype_ranges[] = {152, 48, 34}; // Luma, Chroma DC, Chroma AC
// const int pair_square_ranges[] = {31, 23, 17, 13, 9, 9}; // Luma left, Luma right, ...
const QColor block_colours[] = {Qt::darkBlue, Qt::darkRed, Qt::darkMagenta};

class StegoWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  StegoWidget(QWidget* parent = 0, StegoModel *model = 0);
  
//   QSize sizeHint();
//   QSize minimumSizeHint();
  virtual void updateMinSize() { };
  void updateView();
protected:
  int hd;
  int scale;
  int barHeight;
//   unsigned char **ranges;
  int **ranges;
  StegoModel *model;
public slots:
  void newHorizD(int hd);
  void newScale(int scale);
  void newBarHeight(int bh);
};


class HistogramWidget : public StegoWidget {
//   Q_OBJECT
public:
  HistogramWidget(QWidget *parent = 0, StegoModel *model = 0);
  
  void updateMinSize();
protected:
  void paintEvent(QPaintEvent *event);
  void paintHistogram(int x, int y, double* values, int dim, QPainter *painter);
};


class PairWidget : public StegoWidget {
//   Q_OBJECT
public:
  PairWidget(QWidget *parent = 0, StegoModel *model = 0);
  
  void updateMinSize();
protected:
  QImage *cache;
  
  void paintCache();
  void paintEvent(QPaintEvent *event);
  void paintHistogramBox(QPainter* painter, double* v, int& x, int& y);
  int paintSquares(QPainter* painter, double* v, int& x, int& y, int w, int h);
};

#endif  /*STEGOWIDGETS*/