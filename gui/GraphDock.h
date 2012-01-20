#ifndef STEGOGRAPH
#define STEGOGRAPH

#include <Stegosaurus.h>
#include <QtGui>

class GraphDock : public QDockWidget {
  Q_OBJECT
public:
  GraphDock(QWidget *parent = 0);
};

class StegoTableModel : public QAbstractTableModel, public StegoView {
public:
  StegoTableModel(QObject* parent = 0, StegoModel *stegModel = 0);
    
  int rowCount(const QModelIndex &parent) const;  
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  
  void updateView();
protected:
  int dim;
  double *sigma;
  StegoModel *stegModel;
};

class TableWidget : public QWidget, public StegoView {
  Q_OBJECT
public:
  TableWidget(QWidget* parent = 0, StegoModel *model = 0);
  
  void updateView();
protected:
  QLayout *layout;
  StegoModel *model;
  QTableView *table;
  StegoTableModel *stmodel;
};


class TableDock : public QDockWidget {
  Q_OBJECT
public:
  TableDock(QWidget *parent = 0, StegoModel *model = 0);
protected:
};


#endif  /*STEGOGRAPH*/