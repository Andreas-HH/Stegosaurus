#include "StegoFrame.h"

// extern int openD(StegoModel *model, int num, char **fnames) {
//   for (int i = 0; i < num; i++) {
//     printf("fnames[%i] = %s \n", i, fnames[i]);
//     model->openDirectory(fnames[i]);
//     free(fnames[i]);
//   }
//   free(fnames);
//   model->estimateMus();
// }

StegoFrame::StegoFrame(QWidget* parent) : QMainWindow(parent) {
  QSplitter *central = new QSplitter(Qt::Vertical);
  QScrollArea *hscroll = new QScrollArea();
  QScrollArea *pscroll = new QScrollArea();
  
  model = new StegoModel();
  model->addView(this);
  
  hw = new HistogramWidget(this, model);
  pw = new PairWidget(this, model);

  fdock = new FeatureDock(this, model);
  cdock = new ConfigDock(hw, pw, this);
  gdock = new GraphDock(this);
  tdock = new TableDock(this, model);
  
  hscroll->setWidgetResizable(1);
  hscroll->setAlignment(Qt::AlignLeft);
//   hscroll->setBackgroundRole(QPalette::Dark);
  hscroll->setWidget(hw);
  pscroll->setWidgetResizable(1);
  pscroll->setAlignment(Qt::AlignLeft);
//   pscroll->setBackgroundRole(QPalette::Dark);
  pscroll->setWidget(pw);
  central->addWidget(hscroll);
  central->addWidget(pscroll);
//   scroll1->setWidget(central);
//   central->show();
  
  setCentralWidget(central);
  addDockWidget(Qt::LeftDockWidgetArea, fdock);
  addDockWidget(Qt::RightDockWidgetArea, cdock);
  addDockWidget(Qt::BottomDockWidgetArea, gdock);
  addDockWidget(Qt::BottomDockWidgetArea, tdock);
  
  openAction = new QAction(tr("Open"), this);
  mmdAction  = new QAction(tr("MMDs"), this);
  muAction   = new QAction(tr("Mus"), this);
  
  fdial = new QFileDialog(this);
  fdial->setFileMode(QFileDialog::ExistingFiles);
  fdial->setNameFilter(tr("Features (*.fv)"));
  fdial->setAcceptMode(QFileDialog::AcceptOpen);
  fdial->setDirectory(tr("."));
  
  fileMenu = menuBar()->addMenu(tr("Features"));
  calcMenu = menuBar()->addMenu(tr("Calculate"));
  fileMenu->addAction(openAction);
  calcMenu->addAction(muAction);
  calcMenu->addAction(mmdAction);
  
  progress = new QProgressBar(statusBar());
  statusLabel = new QLabel(tr("Ready."));
  statusBar()->addPermanentWidget(statusLabel, 9);
  statusBar()->addPermanentWidget(progress, 1);
  
  connect(openAction, SIGNAL(triggered(bool)), this, SLOT(openCollection()));
  connect(mmdAction, SIGNAL(triggered(bool)), this, SLOT(calcMMDs()));
  connect(muAction, SIGNAL(triggered(bool)), this, SLOT(calcMus()));
//   connect(cdock, SIGNAL(newHorizD(int)), hw, SLOT(newHorizD(int)));
  
  setMinimumSize(1600, 900);
  setWindowTitle("Stegosaurus");
}

void StegoFrame::openCollection() {
  int i;
  QStringList fileNames;

  if (fdial->exec()) {
    fileNames = fdial->selectedFiles();
    for (i = 0; i < fileNames.size(); i++) {
      model->openFile(fileNames.at(i).toStdString().c_str(), i+1, fileNames.size());
    }
  }
//   model->collectionChanged();
}

void StegoFrame::calcMMDs() {
  model->doMMDs();
}

void StegoFrame::calcMus() {
  model->estimateMus();
}


void StegoFrame::updateProgress(double p) {
//   printf("updating progress %f \n", p);
  progress->setValue((int) (p * 100. + 0.5));
  update();
  hw->update();
  pw->update();
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
