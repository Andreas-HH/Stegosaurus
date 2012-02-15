#include "StegoFrame.h"


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
  
  fdial = new QFileDialog(this);
  fdial->setFileMode(QFileDialog::Directory);
  fdial->setNameFilter(tr("Features (*.fv)"));
  fdial->setAcceptMode(QFileDialog::AcceptOpen);
  fdial->setDirectory(tr("."));
  
  fileMenu = menuBar()->addMenu(tr("Features"));
  fileMenu->addAction(openAction);
  
  progress = new QProgressBar(statusBar());
  statusLabel = new QLabel(tr("Ready."));
  statusBar()->addPermanentWidget(statusLabel, 9);
  statusBar()->addPermanentWidget(progress, 1);
  
  connect(openAction, SIGNAL(triggered()), this, SLOT(openCollection()));
//   connect(cdock, SIGNAL(newHorizD(int)), hw, SLOT(newHorizD(int)));
  
  setMinimumSize(1600, 900);
  setWindowTitle("Stegosaurus");
}

void StegoFrame::openCollection() {
  int i;
  QStringList fileNames;
//   QString fileName = QFileDialog::getOpenFileName(this, tr("Open Features"), ".", tr("Features (*.fv)"));
//   printf("Opening features %s \n", fileName.toStdString().c_str());
//   model->openCollection(fileName.toStdString().c_str());
  if (fdial->exec()) {
    fileNames = fdial->selectedFiles();
    for (i = 0; i < fileNames.size(); i++) {
//       printf("Opening features %s \n", fileNames.at(i).toStdString().c_str());
      model->openDirectory(fileNames.at(i).toStdString().c_str());
    }
  }
  model->estimateMus();
}

void StegoFrame::updateProgress(double p) {
//   printf("updating progress %f \n", p);
  progress->setValue((int) (p * 100.));
  update();
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
