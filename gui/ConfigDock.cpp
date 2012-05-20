#include "ConfigDock.h"


ConfigWidget::ConfigWidget(HistogramWidget *hw, PairWidget *pw, QWidget* parent): QWidget(parent) {
  this->hw = hw;
  this->pw = pw;
  layout = new QGridLayout();
  layout->setVerticalSpacing(0);
  
  hdSlider = new QSlider(Qt::Horizontal);
  hdSlider->setMinimum(1);
  hdSlider->setMaximum(20);
  hdSlider->setTickInterval(1);
  hdSlider->setTickPosition(QSlider::TicksBelow);
  scSlider = new QSlider(Qt::Horizontal);
  scSlider->setMinimum(1);
  scSlider->setMaximum(500);
  scSlider->setTickInterval(10);
  scSlider->setTickPosition(QSlider::TicksBelow);
  bhSlider = new QSlider(Qt::Horizontal);
  bhSlider->setMinimum(1);
  bhSlider->setMaximum(500);
  bhSlider->setTickInterval(10);
  bhSlider->setTickPosition(QSlider::TicksBelow);

  connect(hdSlider, SIGNAL(valueChanged(int)), hw, SLOT(newHorizD(int)));
  connect(scSlider, SIGNAL(valueChanged(int)), hw, SLOT(newScale(int)));
  connect(bhSlider, SIGNAL(valueChanged(int)), hw, SLOT(newBarHeight(int)));
  connect(scSlider, SIGNAL(valueChanged(int)), pw, SLOT(newScale(int)));
  connect(bhSlider, SIGNAL(valueChanged(int)), pw, SLOT(newBarHeight(int)));
  hdSlider->setValue(5);
  scSlider->setValue(50);
  bhSlider->setValue(150);

  layout->addWidget(new QLabel(tr("columns:")), 0, 0);
  layout->addWidget(new QLabel(tr("scale:")), 1, 0);
  layout->addWidget(new QLabel(tr("height:")), 2, 0);
  layout->addWidget(hdSlider, 0, 1);
  layout->addWidget(scSlider, 1, 1);
  layout->addWidget(bhSlider, 2, 1);
  setLayout(layout);
}


ConfigDock::ConfigDock(HistogramWidget *hw, PairWidget *pw, QWidget *parent) : QDockWidget(tr("Settings"),  parent) {
  panel = new ConfigWidget(hw, pw, this);
//   connect(panel, SIGNAL(newHorizD(int)), this, SIGNAL(newHorizD(int)));
  setWidget(panel);
}