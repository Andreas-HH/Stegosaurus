#include "StegoWidgets.h"


StegoWidget::StegoWidget(QWidget* parent, StegoModel* model): QWidget(parent) {
//   setMinimumSize(640, 480);
  this->model = model;
  model->addView(this);
  barHeight = 150;
  scale = 50;
  hd = 5;
}


void StegoWidget::updateView() {
//   setMinimumSize(model->getDimension(), 200);
//   printf("updating view \n");
//   setMinimumSizeHint(2000, 1000);
//   printf("update: %i, %i \n", size().width(), size().height());
//   printf("update: %i, %i \n", sizeHint().width(), sizeHint().height());
//   printf("update: %i, %i \n", minimumSize().width(), minimumSize().height());
//   printf("update: %i, %i \n", minimumSizeHint().width(), minimumSizeHint().height());
  updateMinSize();
  update();
}

void StegoWidget::newHorizD(int hd) {
//   printf("Setting hd %i \n", hd);
  this->hd = hd;
  updateMinSize();
  update();
}

void StegoWidget::newScale(int scale) {
  this->scale = scale;
  updateMinSize();
  update();
}

void StegoWidget::newBarHeight(int bh) {
  this->barHeight = bh;
  updateMinSize();
  update();
}


/*QSize StegoWidget::sizeHint() {
  return QSize(2000, 1000);
}

QSize StegoWidget::minimumSizeHint() {
  return QSize(2000, 1000);
}*/

HistogramWidget::HistogramWidget(QWidget* parent, StegoModel *model): StegoWidget(parent, model) {
}

void HistogramWidget::updateMinSize() {
  if (model->getDimension() != -1) {
    setMinimumSize((model->getDimension()/QP_RANGE + 50)*qMin<int>(hd, QP_RANGE) + 50, (barHeight+50)*((QP_RANGE+hd-1)/hd) + 50);
  }
}

void HistogramWidget::paintEvent(QPaintEvent* event) {
//   printf("painting \n");
  int qx, qy;
  int dim = model->getDimension();
  double *mu = model->getMuVector();
  double *max = model->getMaxVector();
  double *qpHist= model->getQPHist();
  QPainter painter(this);
  QString *string = new QString("%1");
  
//   painter.setBrush(Qt::lightGray);
  painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  
//   printf("painting with dim=%i \n", dim);
  if (dim == -1) return;
  
//   for (i = 1; i < dim; i++)
//     painter.drawLine(i-1, 200-mu[i-1],i, 200-mu[i]);
  for (qx = 0; qx < hd; qx++) {
    for (qy = 0; hd*qy + qx < QP_RANGE; qy++) {
      painter.setPen(Qt::darkCyan);
      if (qpHist != 0) {
        painter.drawText(qx*(dim/QP_RANGE+50)+50+dim/(2*QP_RANGE), (barHeight+50)*(qy+1)+20, string->arg(qpHist[qx+hd*qy]));
      }
      
      if (max != 0) {
        painter.setPen(Qt::gray);
        paintHistogram(qx*(dim/QP_RANGE+50) + 50, (barHeight+50)*(qy+1), &(max[(qx+hd*qy)*dim/QP_RANGE]), dim/QP_RANGE, &painter);
      }
      
      if (mu != 0) {
        painter.setPen(Qt::darkBlue);
        paintHistogram(qx*(dim/QP_RANGE+50) + 50, (barHeight+50)*(qy+1), &(mu[(qx+hd*qy)*dim/QP_RANGE]), dim/QP_RANGE, &painter);
      }
    }
  }
}

void HistogramWidget::paintHistogram(int x, int y, double* values, int dim, QPainter *painter) {
  int i;
  QPen pen = painter->pen();
  
  painter->setPen(Qt::black);
  painter->drawLine(x, y+2, x+dim, y+2);
//   painter->drawLine(x+174,y+3,x+174,y+4);
//   painter->drawLine(x+222,y+3,x+222,y+4);
  painter->drawLine(x+hist_blocktype_ranges[0],y+3,x+hist_blocktype_ranges[0],y+4);
  painter->drawLine(x+hist_blocktype_ranges[0]+hist_blocktype_ranges[1],y+3,x+hist_blocktype_ranges[0]+hist_blocktype_ranges[1],y+4);
  painter->setPen(pen);
  for (i = 0; i < dim; i++) {
    if (values[i] > 0) {   // make zeros easily spottable
//       if (scale*values[i]<barHeight)
// 	painter->drawLine(x+i, y, x+i, y-scale*values[i]);
//       else 
      painter->drawLine(x+i, y, x+i, qMax<int>(y-barHeight, y-scale*values[i]));
    }
  }
/*
  painter.setPen(Qt::darkBlue);
  for (i = 0; i < dim/20; i++) {
    if (max[dim/20*q + i] > 0)
      painter.drawLine(i, barHeight*(q+1),i, barHeight*(q+1)-mu[dim/20*q + i]);
  }*/
}


PairWidget::PairWidget(QWidget *parent, StegoModel *model): StegoWidget(parent, model) {
  hd = 10;
  cache = new QImage(1000, 200, QImage::Format_ARGB32_Premultiplied);
}

void PairWidget::updateMinSize() {
  paintCache();
}

void PairWidget::paintEvent(QPaintEvent *event) {
  printf("painting! \n");
  QPainter painter;
  painter.begin(this);

  painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  painter.drawImage(0, 0, *cache);
  
  painter.end();
}


void PairWidget::paintCache() {
//   printf("painting \n");
  int qx, qy;
  int dim = model->getDimension();
  int **ranges = model->getRanges();
  double *mu = model->getMuVector();
  double *max = model->getMaxVector();
  double *qpHist= model->getQPHist();
  QColor *bg = new QColor(0, 0, 0, 0);
  QPainter painter(cache);
  QString *string = new QString("%1");
  
  if (ranges != 0)
    printf("ranges: %i, %i \n", ranges[0][0], ranges[0][1]);
  
//   painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  cache->fill(bg->rgba());
  if (dim == -1) return;
  
  for (qx = 0; qx < hd; qx++) {
    for (qy = 0; hd*qy + qx < QP_RANGE; qy++) {
//       painter.setPen(Qt::darkCyan);
//       if (qpHist != 0) {
//         painter.drawText(qx*(dim/QP_RANGE+50)+50+dim/(2*QP_RANGE), (barHeight+50)*(qy+1)+20, string->arg(qpHist[qx+hd*qy]));
//       }
      
//       if (max != 0) {
//         painter.setPen(Qt::gray);
//         paintHistogram(qx*(dim/QP_RANGE+50) + 50, (barHeight+50)*(qy+1), &(max[(qx+hd*qy)*dim/QP_RANGE]), dim/QP_RANGE, &painter);
//       }
      
      if (mu != 0) {
        painter.setPen(Qt::darkBlue);
// 	printf("call paint with: %i, %i \n", qx*(pair_square_ranges[1]+pair_square_ranges[3]+pair_square_ranges[5]+40)+20, (pair_square_ranges[0]+20)*(qy+1));
        paintSquares(qx*(pair_square_ranges[1]+pair_square_ranges[3]+pair_square_ranges[5]+40)+20, (pair_square_ranges[0]+20)*(qy)+20, &(mu[(qx+hd*qy)*dim/QP_RANGE]), &painter);
      }
    }
  }
}

void PairWidget::paintSquares(int x, int y, double *values, QPainter *painter) {
  int i, j;
  double *v = values;
  
  painter->setPen(block_colours[0]);
  for (i = 0; i < pair_square_ranges[0]; i++) {
    for (j = 0; j < pair_square_ranges[1]; j++) {
      painter->setOpacity(((double) qMin<int>(barHeight, scale*v[i*pair_square_ranges[1] + j]))/((double) barHeight));
//       printf("drawing point (%i, %i) with opacity %f, scale %i, value %f \n", x+j, y+i, ((double) qMin<int>(barHeight, scale*v[i]))/((double) barHeight), scale);
//       if (i == (pair_square_ranges[0]-1)/2 && j == (pair_square_ranges[1]-1)/2) v--;
      painter->drawPoint(x+j, y+i);
    }
  }
  
  painter->setPen(block_colours[1]);
  v += pair_square_ranges[0]*pair_square_ranges[1];
  for (i = 0; i < pair_square_ranges[2]; i++) {
    for (j = 0; j < pair_square_ranges[3]; j++) {
      painter->setOpacity(((double) qMin<int>(barHeight, scale*v[i*pair_square_ranges[3] + j]))/((double) barHeight));
//       printf("drawing point (%i, %i) with opacity %f, scale %i, value %f \n", x+j, y+i, ((double) qMin<int>(barHeight, scale*v[i]))/((double) barHeight), scale);
//       if (i == (pair_square_ranges[2]-1)/2 && j == (pair_square_ranges[3]-1)/2) v--;
      painter->drawPoint(x+j+28, y+i);
    }
  }
  
  painter->setPen(block_colours[2]);
  v += pair_square_ranges[2]*pair_square_ranges[3];
  for (i = 0; i < pair_square_ranges[4]; i++) {
    for (j = 0; j < pair_square_ranges[5]; j++) {
      painter->setOpacity(((double) qMin<int>(barHeight, scale*v[i*pair_square_ranges[5] + j]))/((double) barHeight));
//       printf("drawing point (%i, %i) with opacity %f, scale %i, value %f \n", x+j, y+i, ((double) qMin<int>(barHeight, scale*v[i]))/((double) barHeight), scale);
//       if (i == (pair_square_ranges[4]-1)/2 && j == (pair_square_ranges[5]-1)/2) v--;
      painter->drawPoint(x+j+32, y+i+21);
    }
  }
  
  painter->setPen(Qt::black);
  painter->setOpacity(0.2);
  painter->drawLine(x-1, y-1, x-1, y+31);
  painter->drawLine(x, y+31, x+22, y+31);
  painter->drawLine(x+23, y+31, x+23, y-1);
  painter->drawLine(x+22, y-1, x, y-1);
  
  painter->drawLine(x+27, y-1, x+27, y+17);
  painter->drawLine(x+28, y+17, x+40, y+17);
  painter->drawLine(x+41, y+17, x+41, y-1);
  painter->drawLine(x+40, y-1, x+28, y-1);
  
  painter->drawLine(x+31, y+20, x+31, y+31);
  painter->drawLine(x+32, y+31, x+40, y+31);
  painter->drawLine(x+41, y+31, x+41, y+20);
  painter->drawLine(x+40, y+20, x+32, y+20);
}

