#ifndef STEGOFDOCK
#define STEGOFDOCK

#include <Stegosaurus.h>
#include <QtGui>


class FeatureWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  FeatureWidget(QWidget* parent = 0, StegoModel *model = 0);
  
  void updateCollection();
protected:
  QTreeView  *fview;
  QStandardItemModel *treeModel;
  QStringList *labels;
  QBoxLayout *layout;
  StegoModel *model;
public slots:
  void selection(QModelIndex mi);
};

class FeatureDock : public QDockWidget {
  Q_OBJECT
public:
  FeatureDock(QWidget *parent = 0, StegoModel *model = 0);
};

#endif  /*STEGOFDOCK*/