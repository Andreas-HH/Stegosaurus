#include "FeatureDock.h"


FeatureWidget::FeatureWidget(QWidget* parent, StegoModel* model): QWidget(parent) {
  this->model = model;
  layout = new QVBoxLayout(this);
  
  model->addView(this);
  fview = new QTreeView(this);
  treeModel = new QStandardItemModel();
  labels = new QStringList();
  *labels << tr("rate [bits/block]") << tr("# vectors") << tr("# files") << tr("dim");
  treeModel->setHorizontalHeaderLabels(*labels);
  fview->setModel(treeModel);

  layout->addWidget(fview);
  setLayout(layout);
}

void FeatureWidget::updateCollection() {
  int i, j;
  FeatureCollection::Iterator *iter;
  featureSet *set;
  QList< QStandardItem* > items;
  QStandardItem *parent;
  
  treeModel->clear();
  treeModel->setHorizontalHeaderLabels(*labels);
  connect(fview, SIGNAL(entered(QModelIndex)), this, SLOT(selection(QModelIndex)));
  connect(fview, SIGNAL(activated(QModelIndex)), this, SLOT(selection(QModelIndex)));
  for (i = 0; i < 10; i++) {
    parent = new QStandardItem(tr("Hans"));
    for (j = 0; j < 8; j++) {
      iter = model->getFeatureIterator(i, j);
      if (iter != 0) {
	while (iter->hasNext()) {
// 	  printf("Have some iterator! \n");
	  set = iter->next();
	  items.append(new QStandardItem(tr("%1").arg(set->rate)));
	  items.append(new QStandardItem(tr("%1").arg(set->M)));
	  items.append(new QStandardItem(tr("%1").arg(set->id)));
	  items.append(new QStandardItem(tr("%1").arg(set->dim)));
	  parent->appendRow(items);
	  items.clear();
	}
	treeModel->appendRow(parent);	
      }
    }
  }
  update();
}

void FeatureWidget::selection(QModelIndex mi) {
  printf("selected: %i, %i \n", mi.column(), mi.row());
  printf("and parent: %i, %i \n", mi.parent().column(), mi.parent().row());
}


FeatureDock::FeatureDock(QWidget *parent, StegoModel *model) : QDockWidget(tr("Collection"),  parent) {
  FeatureWidget *fw = new FeatureWidget(this, model);
  
  setWidget(fw);
}

