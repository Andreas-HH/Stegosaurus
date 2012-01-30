#ifndef STEGOGRAPH
#define STEGOGRAPH

#include <Stegosaurus.h>
#include <QtGui>

class GraphDock : public QDockWidget {
  Q_OBJECT
public:
  GraphDock(QWidget *parent = 0);
};

#endif  /*STEGOGRAPH*/
