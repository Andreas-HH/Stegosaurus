#include "FeatureDock.h"


FeatureWidget::FeatureWidget(QWidget* parent, StegoModel* model): QWidget(parent) {
  layout = new QVBoxLayout(this);
//   progress = new QProgressBar(this);
//   QLabel *label = new QLabel(tr("I am a registered view!"));
//   progress->setMaximum(100);
//   progress->setValue(30);
  
  layout->addStretch();
//   layout->addWidget(progress);
  this->model = model;
  model->addView(this);
  updateView();
  setLayout(layout);
}

/*void FeatureWidget::updateProgress(double p) {
  printf("updating progress %f \n", p);
  progress->setValue((int) (p * 100.));
  update();
}*/


void FeatureWidget::updateCollection() {
//   printf("updating collection \n");
  featureSet *set;
  QLabel *label;
  delete layout;
  layout = new QVBoxLayout(this);
  
//   collection = 0;
  collection = model->getCollection();
  if (collection != 0) {
//     collection->rewind();
    
    while (collection->hasNext()) {
      set = collection->nextFeatureSet();
//       if (collection->isSelected())
// 	label = new QLabel(tr(set->name).append(tr(" (S)")));
//       else 
	label = new QLabel(tr(set->name));
      layout->addWidget(label);
    }
    layout->addStretch();
//     layout->addWidget(progress);
//     collection->rewind();
  }
  update();
}

FeatureDock::FeatureDock(QWidget *parent, StegoModel *model) : QDockWidget(tr("Collection"),  parent) {
  FeatureWidget *fw = new FeatureWidget(this, model);
  
  setWidget(fw);
}

