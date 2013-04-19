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
  int qp_range = model->getQPRange();
  
  if (model->getDimension() != -1) {
    setMinimumSize((model->getDimension()/qp_range + 70)*qMin<int>(hd, qp_range) + 50, (barHeight+50)*((qp_range+hd-1)/hd) + 50);
  }
}

void HistogramWidget::paintEvent(QPaintEvent* event) {
//   printf("painting \n");
  int delta = 0;
  int qx, qy;
  int qp_range = model->getQPRange();
  int dim = model->getDimension();
  int hdim = model->getHistDim();
  int pdim = model->getPairDim();
  int udim = model->getUvsVDim();
  double *mu = model->getMuVector();
  double *max = model->getMaxVector();
//   double *qpHist= model->getQPHist();
  QPainter painter(this);
  QString *string = new QString("%1");
  
//   painter.setBrush(Qt::lightGray);
  painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  
//   printf("painting with dim=%i \n", dim);
  if (dim == -1) return;
  
//   for (i = 1; i < dim; i++)
//     painter.drawLine(i-1, 200-mu[i-1],i, 200-mu[i]);
  for (qx = 0; qx < hd; qx++) {
    for (qy = 0; hd*qy + qx < qp_range; qy++) {
//       painter.setPen(Qt::darkCyan);
//       if (qpHist != 0) {
//         painter.drawText(qx*(dim/qp_range+50)+50+dim/(2*qp_range), (barHeight+50)*(qy+1)+20, string->arg(qpHist[qx+hd*qy]));
//       }
      
      if (max != 0) {
        painter.setPen(Qt::lightGray);
        paintHistogram(qx*(hdim/qp_range+50) + 50, (barHeight+50)*(qy+1), max, hdim/qp_range, &painter);
	paintHistogram(qx*(pdim/qp_range+50) + hdim + 60, (barHeight+50)*(qy+1), max+hdim, pdim/qp_range, &painter);
	paintHistogram(qx*(udim/qp_range+50) + hdim + pdim + 70, (barHeight+50)*(qy+1), max+hdim+pdim, udim/qp_range, &painter);
      }
      
      if (mu != 0) {
        painter.setPen(Qt::darkBlue);
//         paintHistogram(qx*(dim/qp_range+50) + 50, (barHeight+50)*(qy+1), &(mu[(qx+hd*qy)*dim/qp_range]), dim/qp_range, &painter);
// 	paintHistogram(qx*(hdim/qp_range+50) + 50, (barHeight+50)*(qy+1), &(mu[(qx+hd*qy)*hdim/qp_range]), hdim/qp_range, &painter);
// 	paintHistogram(qx*(pdim/qp_range+50) + hdim + 60, (barHeight+50)*(qy+1), &(mu[(qx+hd*qy)*hdim/qp_range]), pdim/qp_range, &painter);
// 	paintHistogram(qx*(udim/qp_range+50) + hdim + pdim + 70, (barHeight+50)*(qy+1), &(mu[(qx+hd*qy)*udim/qp_range]), udim/qp_range, &painter);
	paintHistogram(qx*(hdim/qp_range+50) + 50, (barHeight+50)*(qy+1), mu, hdim/qp_range, &painter);
	paintHistogram(qx*(pdim/qp_range+50) + hdim + 60, (barHeight+50)*(qy+1), mu+hdim, pdim/qp_range, &painter);
	paintHistogram(qx*(udim/qp_range+50) + hdim + pdim + 70, (barHeight+50)*(qy+1), mu+hdim+pdim, udim/qp_range, &painter);
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
//   painter->drawLine(x+hist_blocktype_ranges[0],y+3,x+hist_blocktype_ranges[0],y+4);
//   painter->drawLine(x+hist_blocktype_ranges[0]+hist_blocktype_ranges[1],y+3,x+hist_blocktype_ranges[0]+hist_blocktype_ranges[1],y+4);
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
//   printf("painting! \n");
  QPainter painter;
  painter.begin(this);

  painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  painter.drawImage(0, 0, *cache);
  
  painter.end();
}


void PairWidget::paintCache() {
//   printf("painting \n");
  int i, j;
  int qx, qy;
  int qp_range = model->getQPRange();
  int dim = model->getDimension();
  int hdim = model->getHistDim();
  int pdim = model->getPairDim();
  int udim = model->getUvsVDim();
  int x, y;
  int delta = 0;
//   int **ranges = model->getRanges();
  double *mu = model->getMuVector();
  double *max = model->getMaxVector();
//   double *qpHist= model->getQPHist();
  QColor *bg = new QColor(0, 0, 0, 0);
  QPainter painter(cache);
  QString *string = new QString("%1");
  
  ranges = model->getRanges();
//   printf("qprange: %i \n", qp_range);
//   if (ranges != 0)
//     printf("ranges: %i, %i \n", ranges[0][0], ranges[0][1]);
  
//   painter.fillRect(QRect(0, 0, size().width(), size().height()), Qt::white);
  cache->fill(bg->rgba());
  if (dim == -1 || qp_range == -1) return;
  
  for (qx = 0; qx < hd; qx++) {
    for (qy = 0; hd*qy + qx < qp_range; qy++) {
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
	
	x = 20;
	y = 20;
	for (i = 0; i < qp_range; i++) {
	  paintHistogramBox(&painter, mu + delta, x, y);
	  delta += hdim;
// 	  mu += i*hdim;
	}
	
	for (i = 0; i < qp_range; i++) {
	  painter.setPen(block_colours[0]);
	  delta += paintSquares(&painter, mu + delta, x, y, ranges[0][1], ranges[0][0]);
	  painter.setPen(block_colours[1]);
	  delta += paintSquares(&painter, mu + delta, x, y, ranges[1][1], ranges[1][0]);
	  painter.setPen(block_colours[2]);
	  delta += paintSquares(&painter, mu + delta, x, y, ranges[2][1], ranges[2][0]);
	}
	for (i = 0; i < qp_range; i++) {
	  painter.setPen(Qt::darkYellow);
	  delta += paintSquares(&painter, mu + delta, x, y, ranges[1][0], ranges[1][0]);
	  for (j = 0; j < 15; j++) {
	    if (ranges[2][j] > 1) {
	      delta += paintSquares(&painter, mu + delta, x, y, ranges[2][j], ranges[2][j]);
	    }
	  }
	}
	
// 	printf("delta: %i \n", delta);
	
// 	printf("call paint with: %i, %i \n", qx*(pair_square_ranges[1]+pair_square_ranges[3]+pair_square_ranges[5]+40)+20, (pair_square_ranges[0]+20)*(qy+1));
//         paintSquares(qx*(ranges[1]+pair_square_ranges[3]+pair_square_ranges[5]+40)+20, (pair_square_ranges[0]+20)*(qy)+20, &(mu[(qx+hd*qy)*dim/QP_RANGE]), &painter);
      }
    }
  }
}

void PairWidget::paintHistogramBox(QPainter *painter, double *v, int &x, int &y) {
  int i, j;
  int count = 0;
  int delta;
  
  painter->setPen(block_colours[0]);
  for (i = 0; i < 16; i++) {
    for (j = 0; j < ranges[0][i]; j++) {
      delta = (ranges[0][0] - ranges[0][i])/2;
      painter->setOpacity(sqrt(((double) qMin<int>(barHeight, std::abs(scale*v[count++])))/((double) barHeight)));
      painter->drawPoint(x + j + delta, y + i);
    }
  }
  x += ranges[0][0]+10;
  
  painter->setPen(block_colours[1]);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < ranges[1][i]; j++) {
      delta = (ranges[1][0] - ranges[1][i])/2;
      painter->setOpacity(sqrt(((double) qMin<int>(barHeight, abs(scale*v[count++])))/((double) barHeight)));
      painter->drawPoint(x + j + delta, y + i);
    }
  }
  x += ranges[1][0] + 10;
  
  painter->setPen(block_colours[2]);
  for (i = 0; i < 15; i++) {
    for (j = 0; j < ranges[2][i]; j++) {
      delta = (ranges[2][0] - ranges[2][i])/2;
      painter->setOpacity(sqrt(((double) qMin<int>(barHeight, std::abs(scale*v[count++])))/((double) barHeight)));
      painter->drawPoint(x + j + delta, y + i);
    }
  }
  x += ranges[2][0] + 10;
//   printf("count: %i \n", count);
}

int PairWidget::paintSquares(QPainter *painter, double *v, int &x, int &y, int w, int h) {
  int i, j;
  int count = 0;
  
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      painter->setOpacity(sqrt(((double) qMin<int>(barHeight, std::abs(scale*v[count++])))/((double) barHeight)));
      painter->drawPoint(x+j, y+h-1-i);
    }
  }
  x += ranges[0][1]+10;
  
  return count;
/*  
  painter->setPen(block_colours[1]);
  for (i = 0; i < ranges[1][0]; i++) {
    for (j = 0; j < ranges[1][1]; j++) {
      painter->setOpacity(((double) qMin<int>(barHeight, abs(scale*v[count++])))/((double) barHeight));
//       printf("drawing point (%i, %i) with opacity %f, scale %i, value %f \n", x+j, y+i, ((double) qMin<int>(barHeight, scale*v[i]))/((double) barHeight), scale);
//       if (i == (pair_square_ranges[0]-1)/2 && j == (pair_square_ranges[1]-1)/2) v--;
      painter->drawPoint(x+j, y+ranges[1][0]-1-i);
    }
  }
  x += ranges[1][1]+10;
  
  painter->setPen(block_colours[2]);
  for (i = 0; i < ranges[2][0]; i++) {
    for (j = 0; j < ranges[2][1]; j++) {
      painter->setOpacity(((double) qMin<int>(barHeight, abs(scale*v[count++])))/((double) barHeight));
//       printf("drawing point (%i, %i) with opacity %f, scale %i, value %f \n", x+j, y+i, ((double) qMin<int>(barHeight, scale*v[i]))/((double) barHeight), scale);
//       if (i == (pair_square_ranges[0]-1)/2 && j == (pair_square_ranges[1]-1)/2) v--;
      painter->drawPoint(x+j, y+ranges[2][0]-1-i);
    }
  }
  x += ranges[2][1]+10;*/
//   printf("count: %i \n", count);
//   
//   painter->setPen(Qt::black);
//   painter->setOpacity(0.2);
//   painter->drawLine(x-1, y-1, x-1, y+31);
//   painter->drawLine(x, y+31, x+22, y+31);
//   painter->drawLine(x+23, y+31, x+23, y-1);
//   painter->drawLine(x+22, y-1, x, y-1);
//   
//   painter->drawLine(x+27, y-1, x+27, y+17);
//   painter->drawLine(x+28, y+17, x+40, y+17);
//   painter->drawLine(x+41, y+17, x+41, y-1);
//   painter->drawLine(x+40, y-1, x+28, y-1);
//   
//   painter->drawLine(x+31, y+20, x+31, y+31);
//   painter->drawLine(x+32, y+31, x+40, y+31);
//   painter->drawLine(x+41, y+31, x+41, y+20);
//   painter->drawLine(x+40, y+20, x+32, y+20);
}

