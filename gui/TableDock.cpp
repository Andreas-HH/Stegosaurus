#include "TableDock.h"


StegoTableModel::StegoTableModel(QObject* parent, StegoModel *stegModel): QAbstractTableModel(parent) {
  dim = stegModel->getSigmaDim();
  sigma = stegModel->getSigma();
  diag = stegModel->getDiag();
  if (sigma == 0) dim = -1;
//   dim = -1;
//   sigma = 0;
//   this->stegModel = stegModel;
//   stegModel->addView(this);
}

int StegoTableModel::rowCount(const QModelIndex& parent) const {
  if (dim == -1) return 0;
  return dim;
}

int StegoTableModel::columnCount(const QModelIndex& parent) const{
  if (dim == -1) return 0;
  return dim+1;
}

QVariant StegoTableModel::data(const QModelIndex& index, int role) const {
//   printf("called data \n");
  if (!index.isValid()) return QVariant();
  
  if (role == Qt::TextAlignmentRole) {
    return int(Qt::AlignRight | Qt::AlignVCenter);
  } else if (role == Qt::DisplayRole) {
    if (index.column() == dim) {
//       return QString("%1").arg(diag[index.row()]);
      return QString("Oh-oh!");
    }
    return QString("%1").arg(sigma[index.row() + dim*index.column()]);
  }
  return QVariant();
}

void StegoTableModel::updateView() {
  dim = stegModel->getSigmaDim();
//   printf("sigmadim = %i \n", dim);
  sigma = stegModel->getSigma();
  if (sigma == 0) dim = -1;
}

TableWidget::TableWidget(QWidget* parent, StegoModel *model): QWidget(parent) {
  this->model = model;
  model->addView(this);
  table = new QTableView();
  stmodel = new StegoTableModel(this, model);
  table->setModel(stmodel);
  layout = new QVBoxLayout(this);
  layout->addWidget(table);
}

void TableWidget::updateView() {
  delete stmodel;
  stmodel = new StegoTableModel(this, model);
//   printf("created new model %i, %i \n", model->getDimension(), model->getSigma());
  table->setModel(stmodel);
  update();
}


TableDock::TableDock(QWidget* parent, StegoModel* model): QDockWidget(tr("Table"), parent) {
  TableWidget *tw = new TableWidget(this, model);
  setWidget(tw);
}