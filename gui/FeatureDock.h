#ifndef STEGOFDOCK
#define STEGOFDOCK

#include <Stegosaurus.h>
#include <QtGui>

static const char *blockstrings[8] = {"clean", "L", "C_dc", "LC_dc", "C_ac", "LC_ac", "C", "LC"};
static const char *slice_types[2]  = {"P", "B"};

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
  
    
  featureSet *set1;
  featureSet *set2;
public slots:
  void selection(QModelIndex mi);
};

class FeatureDock : public QDockWidget {
  Q_OBJECT
public:
  FeatureDock(QWidget *parent = 0, StegoModel *model = 0);
};

class FeatureSelectionItem : public QStandardItem {
public:
  FeatureSelectionItem(const QString &text, featureSet *set);
  featureSet *set; // maybe a get method would be nice
protected:
//   StegoModel *model;
//   featureSet *set;
};

#endif  /*STEGOFDOCK*/