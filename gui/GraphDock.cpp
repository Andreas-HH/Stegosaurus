#include "GraphDock.h"


GraphDock::GraphDock(QWidget *parent) : QDockWidget(tr("Graphs"),  parent) {
  QLabel *label = new QLabel("Some graphs");
//   QLayout *layout = new QBoxLayout(QBoxLayout::TopToBottom, this);
  
//   layout->addWidget(label);
//   setLayout(layout);
  setWidget(label);
}


