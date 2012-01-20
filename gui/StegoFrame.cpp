#include "StegoFrame.h"


StegoFrame::StegoFrame(QWidget* parent) : QMainWindow(parent) {
  QSplitter *central = new QSplitter(Qt::Vertical);
//   QScrollArea *scroll1 = new QScrollArea(this);
  
  model = new StegoModel();
  model->addView(this);

  fdock = new FeatureDock(this, model);
  cdock = new ConfigDock(this);
  gdock = new GraphDock(this);
  tdock = new TableDock(this, model);
  
  hw = new HistogramWidget(this, model);
  pw = new PairWidget(this);
  central->addWidget(hw);
  central->addWidget(pw);
//   scroll1->setWidget(central);
//   central->show();
  
  setCentralWidget(central);
  addDockWidget(Qt::LeftDockWidgetArea, fdock);
  addDockWidget(Qt::RightDockWidgetArea, cdock);
  addDockWidget(Qt::BottomDockWidgetArea, gdock);
  addDockWidget(Qt::BottomDockWidgetArea, tdock);
  
  openAction = new QAction(tr("Open"), this);
  
  fileMenu = menuBar()->addMenu(tr("Features"));
  fileMenu->addAction(openAction);
  
  progress = new QProgressBar(statusBar());
  statusLabel = new QLabel(tr("Ready."));
  statusBar()->addPermanentWidget(statusLabel, 9);
  statusBar()->addPermanentWidget(progress, 1);
  
  connect(openAction, SIGNAL(triggered()), this, SLOT(openCollection()));
  
  setWindowTitle("Stegosaurus");
}

void StegoFrame::openCollection() {
//   printf("Opening features \n");
  model->openCollection("data2");
  model->estimateMus();
}

void StegoFrame::updateProgress(double p) {
//   printf("updating progress %f \n", p);
  progress->setValue((int) (p * 100.));
  update();
}

HistogramWidget::HistogramWidget(QWidget* parent, StegoModel *model): QWidget(parent) {
  this->model = model;
  barHeight = 150;
}

void HistogramWidget::paintEvent(QPaintEvent* event) {
//   printf("painting \n");
  int qx, qy;
  int hw = 5;
  int dim = model->getDimension();
  double *mu = model->getMuVector();
  double *max = model->getMaxVector();
  double *qpHist= model->getQPHist();
  QPainter painter(this);
  QString *string = new QString("%1");
  
//   printf("painting with dim=%i \n", dim);
  if (dim == -1) return;
  
//   for (i = 1; i < dim; i++)
//     painter.drawLine(i-1, 200-mu[i-1],i, 200-mu[i]);
  for (qx = 0; qx < hw; qx++) {
    for (qy = 0; hw*qy + qx < QP_RANGE; qy++) {
      if (hw*qy + qx < QP_RANGE) {
	painter.setPen(Qt::darkCyan);
	painter.drawText(qx*(dim/QP_RANGE+50)+100, 200*(qy+1)-180, string->arg(qpHist[qx+4*qy]));
	
	painter.setPen(Qt::gray);
	paintHistogram(qx*(dim/QP_RANGE+50), 200*(qy+1), &(max[(qx+4*qy)*dim/QP_RANGE]), dim/QP_RANGE, &painter);
	
	painter.setPen(Qt::darkBlue);
	paintHistogram(qx*(dim/QP_RANGE+50), 200*(qy+1), &(mu[(qx+4*qy)*dim/QP_RANGE]), dim/QP_RANGE, &painter);
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
  painter->drawLine(x+ranges[0],y+3,x+ranges[0],y+4);
  painter->drawLine(x+ranges[0]+ranges[1],y+3,x+ranges[0]+ranges[1],y+4);
  painter->setPen(pen);
  for (i = 0; i < dim; i++) {
    if (values[i] > 0) {
      if (50*values[i]<150)
	painter->drawLine(x+i, y, x+i, y-50*values[i]);
      else 
	painter->drawLine(x+i, y, x+i, y-150);
    }
  }
/*
  painter.setPen(Qt::darkBlue);
  for (i = 0; i < dim/20; i++) {
    if (max[dim/20*q + i] > 0)
      painter.drawLine(i, barHeight*(q+1),i, barHeight*(q+1)-mu[dim/20*q + i]);
  }*/
}



void HistogramWidget::updateView() {
//   setMinimumSize(model->getDimension(), 200);
//   printf("updating view \n");
  update();
}


PairWidget::PairWidget(QWidget* parent): QWidget(parent) {
}

void PairWidget::updateView() {
  
}


int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
//   void *space = malloc(sizeof(StegoModel)); // new(space)
//   StegoModel *model = new StegoModel();
  StegoFrame *sframe = new StegoFrame();
  
//   sframe->openCollection();
//   model->openCollection("data4");
//   model->estimateMus();
//   model->openCollection("data");
  
  sframe->show();

  return app.exec();
}
