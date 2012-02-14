#include "Stegosaurus.h"


void pathConcat(const char* a, const char* b, char* result) {
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


FeatureCollection::FeatureCollection(featureHeader *h) {
  memcpy(header, h, sizeof(featureHeader));
}

FeatureCollection::~FeatureCollection() {
  int i;
  map< int, featureSet* >::iterator fiter;
  
  for (fiter = collection.begin(); fiter != collection.end(); fiter++) {
    closeFeatureSet(fiter->second);
  }
}

// void FeatureCollection::rewind() {
//   current_set = 0;
// }
// 
// void FeatureCollection::setSelected() {
// //   printf("setting selcted: %i", current_set);
//   selected = current_set;
// }
// 
// int FeatureCollection::isSelected() {
// //   printf("isSelected, %i, %i, %i \n", current_set, selected, current_set == selected);
//   return (current_set == selected);
// }

int FeatureCollection::addFeatureFile(const char *path, featureHeader *header) {
  int bin = (int) (header->rate/BIN_WIDTH + 0.5);
  featureSet *set;

  printf("adding feature file \n");
  if (collection[bin] == 0) {
    printf("no collection for this yet! \n");
    set = openFeatureSet(path);
    collection[bin] = set;
  } else {
    newFeatureFile(collection[bin], path);
  }
}


int FeatureCollection::getNumSets() {
   return collection.size();
}

int FeatureCollection::getCurrentSet() {
  return current_set;
}

featureSet* FeatureCollection::getFeatureSet(int index) {
//   return collection[index];
  return 0;
}

// // void FeatureCollection::setSelected(int index) {
// //   selected = index;
// // }


featureSet* FeatureCollection::nextFeatureSet() {
//   return collection[current_set++];
  return 0;
}

int FeatureCollection::hasNext() {
  if (current_set < num_sets) return 1;
  else return 0;
}


StegoModel::StegoModel() {
  int i;
//   features = new FeatureCollection();
  current_view = 0;
  steg = init_stego();
  for (i = 0; i < 10; i++) {
    collections[i] = 0;
  }
//   fcol = 0;
//   fcol = new FeatureCollection("");
}

StegoModel::~StegoModel() {
//   printf("destroying model \n");
  close_stego(steg);
}

void StegoModel::estimateMus() {
//   int i;
//   
// //   while (fcol->hasNext()) {
// //   printf("estimate mus, %i \n", fcol->getNumSets());
//   for (i = 0; i < fcol->getNumSets(); i++) {
// //     steg->features = fcol->nextFeatureSet();
//     steg->features = fcol->getFeatureSet(i);
// //     printf("set steg->features \n");
// //     fcol->setSelected(i);
// //     printf("set selected correctly \n");
//     estimateMu(steg);
//     progressChanged(((double) i)/(double) fcol->getNumSets());
//     progressChanged(((double) fcol->getCurrentSet())/(double) fcol->getNumSets());
//   }
// //   estimateSigma(steg);
// //   progressChanged(0.99);
// //   qrHouseholder(steg);
// //   qrUnitTest();
// //   fcol->rewind();
//   modelChanged();
//   progressChanged(0.);
// //   printf("done. \n");
}


void StegoModel::addView(StegoView *view) {
//   views[current_view] = view;
  views.push_back(view);
  current_view++;
}

void StegoModel::modelChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
//   for (i = 0; i < current_view; i++) {
//     views[i]->updateView();
//   }
  for (siter = views.begin(); siter != views.end(); siter++) {
    printf("updateing some view \n");
    (*siter)->updateView();
  }
}

void StegoModel::collectionChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
//   for (i = 0; i < current_view; i++) {
//     views[i]->updateCollection();
//   }
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateCollection();
  }
}

void StegoModel::progressChanged(double p) {
  int i;
  list< StegoView* >::iterator siter;
  
//   for (i = 0; i < current_view; i++) {
//     views[i]->updateProgress(p);
//   }
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateProgress(p);
  }
}

void StegoModel::openDirectory(const char* path) {
  int i;
  int num_sets = 0;
  int bin;
  char *str = (char*) malloc(512*sizeof(char));
  DIR *root = opendir(path);
  FILE *file;
  featureHeader header;
  struct dirent *entry;
  
  if (root == NULL) {
//     printf("root is NULL \n");
    return;
  }
  printf("Root not NULL. \n");
  while ((entry = readdir(root)) != NULL) {
    if (strstr(entry->d_name, ".fv") != NULL) {
      num_sets++;
    }
  }
  rewinddir(root);
  for(i = 0; i < num_sets; ) {
    entry = readdir(root);
    if (strstr(entry->d_name, ".fv") != NULL) {
      pathConcat(path, entry->d_name, str);
      printf("about to open str (=%s) \n", str);
      file = fopen(str, "r");
      printf("about to read header \n");
      readHeader(file, &header);
      printf("successfully read header \n");
      fclose(file);
//       set->name = (char*) malloc(strlen(entry->d_name)*sizeof(char));
//       strcpy(set->name, entry->d_name);
      if (collections[header.method] == 0)
	collections[header.method] = new FeatureCollection(&header);
      collections[header.method]->addFeatureFile(str, &header);
      i++;
    }
  }
//   printf("%i sets, collection[0]->M = %i \n", num_sets, collection[0]->M);
  closedir(root);
  free(str);
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
//   return fcol; // fcol
  return 0;
}

double* StegoModel::getSigma() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->qr;
}

double* StegoModel::getDiag() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->qr_diag;
}

