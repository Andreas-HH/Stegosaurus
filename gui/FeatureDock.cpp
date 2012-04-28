#include "FeatureDock.h"


FeatureWidget::FeatureWidget(QWidget* parent, StegoModel* model): QWidget(parent) {
  this->model = model;
  layout = new QVBoxLayout(this);
  
  model->addView(this);
  fview = new QTreeView(this);
  treeModel = new QStandardItemModel();
  labels = new QStringList();
  *labels << tr("p_hide") << tr("# vectors") << tr("dimension") << tr("# files") << tr("type") << tr("qp offset") << tr("qp range");
  treeModel->setHorizontalHeaderLabels(*labels);
  fview->setModel(treeModel);

  layout->addWidget(fview);
  setLayout(layout);
  
  connect(fview, SIGNAL(entered(QModelIndex)), this, SLOT(selection(QModelIndex)));
  connect(fview, SIGNAL(activated(QModelIndex)), this, SLOT(selection(QModelIndex)));
  
  set1 = 0;
  set2 = 0;
}

void FeatureWidget::updateCollection() {
  int i, j, idx = 0;
  int currentMethod = -1; // some illegal value
//   bool appendM = false;
  FeatureCollection::Iterator *iter;
  featureSet *set;
  QList< QStandardItem* > items;
  FeatureSelectionItem *parentM, *parentA;
  StegoModel::Iterator *miter;
  
  treeModel->clear();
  treeModel->setHorizontalHeaderLabels(*labels);

  if ((set = model->getCleanSet()) != 0) {
//     set1 = set;
    items.append(new FeatureSelectionItem(tr("%1").arg(0), set));
    items.append(new FeatureSelectionItem(tr("%1").arg(set->M), set));
    items.append(new FeatureSelectionItem(tr("%1").arg(set->dim), set));
    items.append(new FeatureSelectionItem(tr("%1").arg(set->num_files), set));
    items.append(new FeatureSelectionItem(tr("%1").arg(slice_types[(int)set->header->slice_type]), set));
    items.append(new FeatureSelectionItem(tr("%1").arg((int) set->header->qp_offset), set));
    items.append(new FeatureSelectionItem(tr("%1").arg((int) set->header->qp_range), set));
    treeModel->appendRow(items);
    items.clear();
  }
  
  miter = model->iterator();
  while(miter->hasNext()) {
    iter = miter->next()->iterator();
//     printf("adding some stego set \n");
    if (miter->getX().first != currentMethod) {
      if (currentMethod != -1)
	treeModel->appendRow(parentM);
      currentMethod = miter->getX().first;
      parentM = new FeatureSelectionItem(tr("%1").arg(miter->getX().first), 0);
    }
    parentA = new FeatureSelectionItem(tr("%1").arg(blockstrings[miter->getX().second]), 0);
    while (iter->hasNext()) {
      set = iter->next();
      items.append(new FeatureSelectionItem(tr("%1").arg(set->prob), set));
      items.append(new FeatureSelectionItem(tr("%1").arg(set->M), set));
      items.append(new FeatureSelectionItem(tr("%1").arg(set->dim), set));
      items.append(new FeatureSelectionItem(tr("%1").arg(set->num_files), set));
      items.append(new FeatureSelectionItem(tr("%1").arg(slice_types[(int)set->header->slice_type]), set));
      items.append(new FeatureSelectionItem(tr("%1").arg((int) set->header->qp_offset), set));
      items.append(new FeatureSelectionItem(tr("%1").arg((int) set->header->qp_range), set));
      parentA->appendRow(items);
      items.clear();
    }
    parentM->appendRow(parentA);
  }
  if (currentMethod != -1)
    treeModel->appendRow(parentM);
  update();
}

void FeatureWidget::selection(QModelIndex mi) {
  FeatureSelectionItem *select = (FeatureSelectionItem*) treeModel->itemFromIndex(mi);
  
//   printf("picking set with p_hide = %g \n", select->set->prob);
  if (select->set != 0)
    model->setFeatures(select->set);
//   printf("so: %s \n", treeModel->itemFromIndex(mi)->text().toStdString().c_str());
//   printf("selected: %i, %i \n", mi.column(), mi.row());
//   printf("and parent: %i, %i \n", mi.parent().column(), mi.parent().row());
}

FeatureDock::FeatureDock(QWidget *parent, StegoModel *model) : QDockWidget(tr("Collection"),  parent) {
  FeatureWidget *fw = new FeatureWidget(this, model);
  
  setWidget(fw);
}

FeatureSelectionItem::FeatureSelectionItem(const QString& text, featureSet* set) : QStandardItem(text) {
  this->set = set;
}


