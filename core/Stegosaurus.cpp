#include "Stegosaurus.h"


void pathConcat(char* a, char* b, char *result) {
  int i;
  int count;
  
  for (i = 0; a[i] != '\0'; i++)
    result[i] = a[i];
  result[i] = '/';
  count = i+1;
  for (i = 0; b[i] != '\0'; i++)
    result[count+i] = b[i];
  result[count+i] = '\0';
  
//   printf(" __ %s __ \n", result);
//   return result;
}


FeatureCollection::FeatureCollection(char* path) {
  int i;
//   int numsets = 0;
  char *str = (char*) malloc(512*sizeof(char));
  
  selected = -1;
  current_set = 0;
  DIR *root = opendir(path);
  struct dirent *entry;
  
  num_sets = 0;
  if (root == NULL) {
//     printf("root is NULL \n");
    return;
  }
//   printf("Root not NULL. \n");
  while ((entry = readdir(root)) != NULL) {
    if (strstr(entry->d_name, ".fv") != NULL) {
      num_sets++;
    }
  }
  collection = (featureSet**) malloc(num_sets*sizeof(featureSet*));
  rewinddir(root);
  for(i = 0; i < num_sets; ) {
    entry = readdir(root);
//     printf("%s \n", pathConcat(path, entry->d_name));
    if (strstr(entry->d_name, ".fv") != NULL) {
//       printf("%s \n", entry->d_name);
      pathConcat(path, entry->d_name, str);
//       printf(" %s \n", str);
      collection[i] = openFeatureSet(str);
      collection[i]->name = (char*) malloc(strlen(entry->d_name)*sizeof(char));
      strcpy(collection[i]->name, entry->d_name);
      i++;
    }
  }
//   printf("%i sets, collection[0]->M = %i \n", num_sets, collection[0]->M);
  closedir(root);
  rewind();
  free(str);
//   printf("FeatureCollection: done. \n");
}

FeatureCollection::~FeatureCollection() {
  int i;
  
//   printf("Closing FeatureCollection. \n");
  for (i = 0; i < num_sets; i++) {
    closeFeatureSet(collection[i]);
  }
  free(collection);
}

void FeatureCollection::rewind() {
  current_set = 0;
}

void FeatureCollection::setSelected() {
//   printf("setting selcted: %i", current_set);
  selected = current_set;
}

int FeatureCollection::isSelected() {
//   printf("isSelected, %i, %i, %i \n", current_set, selected, current_set == selected);
  return (current_set == selected);
}

int FeatureCollection::getNumSets() {
  return num_sets;
}

int FeatureCollection::getCurrentSet() {
  return current_set;
}


featureSet* FeatureCollection::nextFeatureSet() {
  return collection[current_set++];
}

int FeatureCollection::hasNext() {
  if (current_set < num_sets) return 1;
  else return 0;
}


StegoModel::StegoModel() {
//   features = new FeatureCollection();
  current_view = 0;
  steg = init_stego();
  fcol = 0;
//   fcol = new FeatureCollection("");
}

StegoModel::~StegoModel() {
//   printf("destroying model \n");
  close_stego(steg);
}

void StegoModel::estimateMus() {
  while (fcol->hasNext()) {
    steg->features = fcol->nextFeatureSet();
    fcol->setSelected();
    estimateMu(steg);
    progressChanged(((double) fcol->getCurrentSet())/(double) fcol->getNumSets());
  }
  estimateSigma(steg);
  fcol->rewind();
  modelChanged();
  progressChanged(0.);
}


void StegoModel::addView(StegoView *view) {
  views[current_view] = view;
  current_view++;
}

void StegoModel::modelChanged() {
  int i;
  
  for (i = 0; i < current_view; i++) {
    views[i]->updateView();
  }
}

void StegoModel::collectionChanged() {
  int i;
  
  for (i = 0; i < current_view; i++) {
    views[i]->updateCollection();
  }
}

void StegoModel::progressChanged(double p) {
  int i;
  
  for (i = 0; i < current_view; i++) {
    views[i]->updateProgress(p);
  }
}

void StegoModel::openCollection(char* path) {
  if (fcol != NULL) delete fcol;
//   printf("opening collection \n");
//   if (fcol == 0) printf("features is NULL \n");
//   else printf("features not NULL \n");
  fcol = new FeatureCollection(path);
//   if (features->hasNext())
//     steg->features = features->collection[0];
//   modelChanged();
  collectionChanged();
  modelChanged();
//   printf("opened collection! \n");
}

int StegoModel::getDimension() {
  if (steg->features == NULL)
    return -1;
  return steg->features->dim;
}


double* StegoModel::getMaxVector() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->max_vec;
}

double* StegoModel::getMuVector() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->mu;
}

double* StegoModel::getQPHist() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->qp_vec;
}

FeatureCollection* StegoModel::getCollection() {
  return fcol; // fcol
}

double* StegoModel::getSigma() {
    if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->sigma;
}


