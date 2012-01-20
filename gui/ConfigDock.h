#ifndef STEGOCONF
#define STEGOCONF

#include <Stegosaurus.h>
#include <QtGui>


class ConfigDock : public QDockWidget {
  Q_OBJECT
public:
  ConfigDock(QWidget *parent = 0);
};

#endif  /*STEGOCONF*/