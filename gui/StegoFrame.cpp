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

LoadDialog::LoadDialog(QWidget* parent, Qt::WindowFlags f) : QDialog(parent, f) {
  int i;
  
  QVBoxLayout *mainLayout = new QVBoxLayout();
  QHBoxLayout *boxLayout = new QHBoxLayout();
  QHBoxLayout *buttonLayout = new QHBoxLayout();
  
  minMethod = new QComboBox(this);
  maxMethod = new QComboBox(this);
  qpOffset = new QComboBox(this);
  qpRange = new QComboBox(this);
  type = new QComboBox(this);
  ok = new QPushButton(tr("OK"));
  cancel = new QPushButton(tr("Cancel"));
  
  for (i = 1; i < 30; i++) {
    minMethod->addItem(QString("%1").arg(i));
    maxMethod->addItem(QString("%1").arg(i));
  }
  qpOffset->addItem(QString("16"));
  qpOffset->addItem(QString("20"));
  qpOffset->addItem(QString("24"));
  qpOffset->addItem(QString("28"));
  qpRange->addItem(QString("1"));
  type->addItem("P");
  type->addItem("B");
  
  boxLayout->addWidget(minMethod);
  boxLayout->addWidget(maxMethod);
  boxLayout->addWidget(qpOffset);
  boxLayout->addWidget(qpRange);
  boxLayout->addWidget(type);
  
  buttonLayout->addStretch();
  buttonLayout->addWidget(ok);
  buttonLayout->addWidget(cancel);
  buttonLayout->addStretch();
  
  mainLayout->addLayout(boxLayout);
  mainLayout->addStretch();
  mainLayout->addLayout(buttonLayout);
  
//   connect(ok, SIGNAL(pressed()), this, SIGNAL(accepted()));
  connect(ok, SIGNAL(pressed()), this, SLOT(accept()));
  connect(cancel, SIGNAL(pressed()), this, SLOT(reject()));
  
  setLayout(mainLayout);
  setWindowTitle(tr("Load Features"));
}

QString LoadDialog::getMinMethod() {
  return minMethod->currentText();
}

QString LoadDialog::getMaxMethod() {
  return maxMethod->currentText();
}

QString LoadDialog::getQPOffset() {
  return qpOffset->currentText();
}

QString LoadDialog::getQPRange() {
  return qpRange->currentText();
}

QString LoadDialog::getType() {
  return type->currentText();
}

StegoFrame::StegoFrame(QWidget* parent) : QMainWindow(parent) {
  QSplitter *central = new QSplitter(Qt::Vertical);
  QScrollArea *hscroll = new QScrollArea();
  QScrollArea *pscroll = new QScrollArea();
  QDomElement stegoDom;
  
  model = new StegoModel();
  model->addView(this);
  
  hw = new HistogramWidget(this, model);
  pw = new PairWidget(this, model);

  fdock = new FeatureDock(this, model);
  cdock = new ConfigDock(hw, pw, this);
//   gdock = new GraphDock(this);
//   tdock = new TableDock(this, model);
  ldial = new LoadDialog(this);
  
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
  addDockWidget(Qt::BottomDockWidgetArea, cdock);
//   addDockWidget(Qt::BottomDockWidgetArea, gdock);
//   addDockWidget(Qt::BottomDockWidgetArea, tdock);
  
  openFeaturesAction = new QAction(tr("Open"), this);
  loadFeaturesAction = new QAction(tr("Load"), this);
  openDocAction = new QAction(tr("Open"), this);
  saveDocAction = new QAction(tr("Save"), this);
  mmdAction  = new QAction(tr("MMDs"), this);
  svmAction  = new QAction(tr("SVM classification"), this);
  muAction   = new QAction(tr("Mus"), this);
  
  fdial = new QFileDialog(this);
  fdial->setFileMode(QFileDialog::ExistingFiles);
  fdial->setNameFilter(tr("Features (*.fv)"));
  fdial->setAcceptMode(QFileDialog::AcceptOpen);
  fdial->setDirectory(tr("."));
  
  docMenu = menuBar()->addMenu(tr("Document"));
  fileMenu = menuBar()->addMenu(tr("Features"));
  calcMenu = menuBar()->addMenu(tr("Calculate"));
  fileMenu->addAction(openFeaturesAction);
  fileMenu->addAction(loadFeaturesAction);
  docMenu->addAction(openDocAction);
  docMenu->addAction(saveDocAction);
  calcMenu->addAction(muAction);
  calcMenu->addAction(mmdAction);
  calcMenu->addAction(svmAction);
  
  progress = new QProgressBar(statusBar());
  statusLabel = new QLabel(tr("Ready."));
  statusBar()->addPermanentWidget(statusLabel, 9);
  statusBar()->addPermanentWidget(progress, 1);
  
  connect(openFeaturesAction, SIGNAL(triggered(bool)), this, SLOT(openCollection()));
  connect(loadFeaturesAction, SIGNAL(triggered(bool)), ldial, SLOT(open()));
  connect(saveDocAction, SIGNAL(triggered(bool)), this, SLOT(saveXML()));
  connect(mmdAction, SIGNAL(triggered(bool)), this, SLOT(calcMMDs()));
  connect(svmAction, SIGNAL(triggered(bool)), this, SLOT(classify()));
  connect(muAction, SIGNAL(triggered(bool)), this, SLOT(calcMus()));
  connect(ldial, SIGNAL(accepted()), this, SLOT(loadFeatures()));
  
  document = new QDomDocument();
  xmlFile = new QFile("stegodoc.xml");
  if (xmlFile->exists()) {
    openXML();
//     if (document->firstChildElement().tagName() != QString("stegosaurus"))
//       printf("lol \n");
    sets = document->firstChildElement().firstChildElement(QString("sets"));
    if (sets.isNull()) {
      sets = document->createElement(QString("sets"));
      stegoDom.appendChild(sets);
    }
    mmds = document->firstChildElement().firstChildElement(QString("mmds"));
    if (sets.isNull()) {
      sets = document->createElement(QString("mmds"));
      stegoDom.appendChild(sets);
    }
  } else {
    stegoDom = document->createElement(QString("stegosaurus"));
    sets = document->createElement(QString("sets"));
    stegoDom.appendChild(sets);
    document->appendChild(stegoDom);
  }
  
  setMinimumSize(1280, 720);
  setWindowTitle("Stegosaurus");
}

void StegoFrame::openCollection() {
  int i;
  QStringList fileNames;
  featureHeader header;
  QDomElement currentSet;// = sets.firstChildElement(QString("featureSet"));
  QDomElement currentFile;
  QDomText text;

  if (fdial->exec()) {
    fileNames = fdial->selectedFiles();
    for (i = 0; i < fileNames.size(); i++) {
      if (model->openFile(fileNames.at(i).toStdString().c_str(), i+1, fileNames.size(), header) == 0) {
        for (currentSet = sets.firstChildElement(QString("featureSet")); !currentSet.isNull(); currentSet = currentSet.nextSiblingElement(QString("featureSet"))) {
    // 	  printf("checking method %i \n", currentSet.attribute(tr("method")).toInt());
          if (currentSet.attribute(tr("method")).toInt() == (int) header.method &&
                currentSet.attribute(tr("qp_range")).toInt() == (int) header.qp_range &&
                currentSet.attribute(tr("qp_offset")).toInt() == (int) header.qp_offset &&
                currentSet.attribute(tr("type")) == slice_types[(int)header.slice_type]) {
    // 	    printf("found some set in document \n");
            break;
          }
        }
        if (currentSet.isNull()) { // need to insert new featureSet
              currentSet = document->createElement(QString("featureSet"));
          currentSet.setAttribute(QString("method"), QString("%1").arg((int)header.method));
          currentSet.setAttribute(QString("type"), QString(slice_types[(int)header.slice_type]));
    // 	  currentSet.setAttribute(QString("dimension"), QString("%1").arg((int)header.dim));
          currentSet.setAttribute(QString("qp_offset"), QString("%1").arg((int)header.qp_offset));
          currentSet.setAttribute(QString("qp_range"), QString("%1").arg((int)header.qp_range));
          sets.appendChild(currentSet);
    // 	  currentSet = currentSet;
        }
        for (currentFile = currentSet.firstChildElement(QString("file")); !currentFile.isNull(); currentFile = currentFile.nextSiblingElement(QString("file"))) {
    // 	  printf("currently looking at: %s \n", currentFile.text().toStdString().c_str());
          if (currentFile.text().compare(fileNames.at(i)) == 0)
            break;
          // compare here to avoid duplicates
        }
        if (currentFile.isNull()) {
    // 	  printf("NOT FOUND! (%s) \n", fileNames.at(i).toStdString().c_str());
          currentFile = document->createElement(QString("file"));
          text = document->createTextNode(fileNames.at(i));
          currentFile.appendChild(text);
          currentSet.appendChild(currentFile);
        }
      }
    }
  }
  model->collectionChanged();
}

void StegoFrame::saveXML() {
//   QFile *output = new QFile("output.xml");
  xmlFile->open(QFile::WriteOnly | QFile::Text);
  QTextStream *out = new QTextStream(xmlFile);
  
  document->save(*out, 2);
  xmlFile->close();
}

void StegoFrame::openXML() {
  int errorLine, errorColumn;
  QString errorStr;
//   QFile file("output.xml");
  
  xmlFile->open(QFile::ReadOnly | QFile::Text);
  document->setContent(xmlFile, false, &errorStr, &errorLine, &errorColumn);
  xmlFile->close();
}

void StegoFrame::loadFeatures() {
  QDomElement currentSet;
  QDomElement currentFile;
  featureHeader header;
  int minMethod = ldial->getMinMethod().toInt();
  int maxMethod = ldial->getMaxMethod().toInt();
  int currentMethod;
  QList<QString> fileList;
  QString f;
  
  for (currentSet = sets.firstChildElement(QString("featureSet")); !currentSet.isNull(); currentSet = currentSet.nextSiblingElement(QString("featureSet"))) {
    currentMethod = currentSet.attribute(QString("method")).toInt();
    if (currentMethod == 0 || (currentMethod >= minMethod && currentMethod <= maxMethod)) {
      if (currentSet.attribute(QString("qp_offset")) == ldial->getQPOffset() && 
          currentSet.attribute(QString("qp_range")) == ldial->getQPRange() &&
          currentSet.attribute(QString("type")) == ldial->getType()) {
        for (currentFile = currentSet.firstChildElement(QString("file")); !currentFile.isNull(); currentFile = currentFile.nextSiblingElement(QString("file"))) {
//           printf("Loading something: %s \n", currentFile.text().toStdString().c_str());
          fileList.append(currentFile.text());
//           model->openFile(currentFile.text().toStdString().c_str(), 0, 1, header);
        }
      }
    }
  }
  qSort(fileList);
  foreach (f, fileList) {
    printf("Loading something: %s \n", f.toStdString().c_str());
    model->openFile(f.toStdString().c_str(), 0, 1, header);
  }
  model->collectionChanged();
}

void StegoFrame::calcMMDs() {
  model->doMMDs();
  
}

void StegoFrame::calcMus() {
  model->estimateMus();
}

void StegoFrame::classify() {
  model->runClassifier();
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
