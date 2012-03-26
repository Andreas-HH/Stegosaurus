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
}


FeatureCollection::FeatureCollection(featureHeader *h) {
  memcpy(&header, h, sizeof(featureHeader));
//   cleanSet = 0;
}

FeatureCollection::~FeatureCollection() {
  int i;
  map< int, featureSet* >::iterator fiter;
  
  for (fiter = collection.begin(); fiter != collection.end(); fiter++) {
    closeFeatureSet(fiter->second);
  }
}

int FeatureCollection::addFeatureFile(const char* path, featureHeader* header, stegoContext* steg) {
  int bin;// = (int) (header->rate/BIN_WIDTH + 0.5);
//   printf("rate = %g, bin = %i \n", header->rate, bin);
  featureSet *set;

  // there should be some health-checking here!
  
  printf("adding feature file \n");
//   if (header->method == 0) {
//     if (cleanSet == 0) {
//       cleanSet = openFeatureSet(path);
//     } else {
//       newFeatureFile(cleanSet, path);
//     }
//   } else {
    bin = (int) (header->rate/BIN_WIDTH + 0.5);
    if (collection[bin] == 0) {
//       printf("no collection for this yet! \n");
      set = openFeatureSet(path, steg);
      set->rate = header->rate;
      set->id = bin;
      collection[bin] = set;
//       printf("just defined collection[%i] \n", bin);
    } else {
//       printf("There already is something \n");
  //     newFeatureFile(collection[bin], path);
    }
//   }
//   printf("added new feature file \n");
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

FeatureCollection::Iterator* FeatureCollection::iterator() {
  return new FeatureCollection::Iterator(this);
}


FeatureCollection::Iterator::Iterator(FeatureCollection* f) {
  fc = f;
  iter = f->collection.begin();
}

bool FeatureCollection::Iterator::hasNext() {
  return (iter != fc->collection.end());
}

featureSet* FeatureCollection::Iterator::next() {
  featureSet *set = iter->second;
  iter++;
  return set;
}


StegoModel::StegoModel() {
  int i, j;
//   features = new FeatureCollection();
//   current_view = 0;
  ranges = 0;
  steg = init_stego();
  cleanSet = 0;
  mc = 0;
//   for (i = 0; i < 10; i++) {
//     for (j = 0; j < 8; j++) {
//       collections[i][j] = 0;
//     }
//   }
}

StegoModel::~StegoModel() {
  if (mc != 0) closeMMD(*mc);
  close_stego(steg);
}

void StegoModel::estimateMus() {
  int i, j;
  FeatureCollection::Iterator* iter;
  
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 8; j++) {
      if (collections[i][j] != 0) {
	iter = collections[i][j]->iterator();
	while (iter->hasNext()) {
	  steg->features = iter->next();
	  printf("About to estimate mu with dim=%i, M=%i \n", steg->features->dim, steg->features->M);
	  estimateMu(steg);
	}
      }
    }
  }
  modelChanged();
}

// we don't want to give sets directly, rather indices inside the collection
double StegoModel::doMMD(featureSet *clean, featureSet *stego) {
//   mmdContext mc;
  if (mc == 0) {
    mc = (mmdContext*) malloc(sizeof(mmdContext));
    mc->clean = clean;
    mc->stego = stego;
    initMMD(steg, *mc);
    estimateGamma(steg, *mc);
  }
  mc->stego = stego;
//   mc->stego = mc->clean;
  estimateMMD(steg, *mc);
//   printf("used gamma: %g \n", mc->gamma);
  
  return mc->mmd;
}


void StegoModel::setFeatures(featureSet* set) {
  steg->features = set;
  modelChanged();
}


void StegoModel::addView(StegoView *view) {
  views.push_back(view);
//   current_view++;
}

void StegoModel::modelChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
    printf("updateing some view \n");
    (*siter)->updateView();
  }
}

void StegoModel::collectionChanged() {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateCollection();
  }
}

void StegoModel::progressChanged(double p) {
  int i;
  list< StegoView* >::iterator siter;
  
  for (siter = views.begin(); siter != views.end(); siter++) {
    (*siter)->updateProgress(p);
  }
}

void StegoModel::openDirectory(const char* path) {
  int i, j, k;
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
//   printf("Root not NULL. \n");
  while ((entry = readdir(root)) != NULL) {
    if (strstr(entry->d_name, ".fv") != NULL) {
      num_sets++;
    }
  }
//   printf("about to open dir with %i feature files \n", num_sets);
  rewinddir(root);
  for(i = 0; i < num_sets; ) {
    entry = readdir(root);
    if (strstr(entry->d_name, ".fv") != NULL) {
      pathConcat(path, entry->d_name, str);
//       printf("about to open str (=%s) \n", str);
      file = fopen(str, "r");
//       printf("about to read header \n");
      readHeader(file, &header);
      if (ranges == 0) {
	printf("first file is being added! \n");
	ranges = (int**) malloc(3*sizeof(int*));
	if (header.pair) {
	  for (k = 0; k < 3; k++) {
	    ranges[k] = (int*) malloc((2*header.ranges[k][0]+1)*sizeof(int));
	    for (j = 0; j < 2*header.ranges[k][0]+1; j++)
	      ranges[k][j] = header.ranges[k][1];
	  }
	} else {
	  for (k = 0; k < 3; k++) {
	    ranges[k] = (int*) malloc(num_coefs[k]*sizeof(int));
	    for (j = 0; j < num_coefs[k]; j++)
	      ranges[k][j] = header.ranges[k][j];
	  }
	}
      }
//       printf("successfully read header \n");
      fclose(file);
//       set->name = (char*) malloc(strlen(entry->d_name)*sizeof(char));
//       strcpy(set->name, entry->d_name);
      printf("method: %i \n", header.method);
      printf("qp range: %i \n", header.qp_range);
      if (header.method == 0) {
	if (cleanSet == 0) {
	  printf("cleanSet shouldn't be null from now on! %s ", str);
	  cleanSet = openFeatureSet(str, steg);
	  if (cleanSet != 0)
	    printf("indeed! \n");
	} else {
	  newFeatureFile(cleanSet, str);
	}
      } else {
	if (collections[header.method][header.accept] == 0) {
	  collections[header.method][header.accept] = new FeatureCollection(&header);
	  printf("Created new collection for method %i and accept %i \n", header.method, header.accept);
	}
        collections[header.method][header.accept]->addFeatureFile(str, &header, steg);
      }
/*      printf("[%i][%i][%i][%i] = %i \n", header.video_bitrate, header.pair, header.method, header.accept, collections[header.video_bitrate][header.pair][header.method][header.accept]);
      if ((((collections[header.video_bitrate])[header.pair])[header.method])[header.accept] == 0) {
	collections[header.video_bitrate][header.pair][header.method][header.accept] = new FeatureCollection(&header);
	printf("Created new collection for method %i and accept %i \n", header.method, header.accept);
      }
      collections[header.video_bitrate][header.pair][header.method][header.accept]->addFeatureFile(str, &header);*/
      i++;
    }
  }
  collectionChanged();
  
//   // print results
//   FeatureCollection::Iterator* iter;
//   printf("printing collections: \n");
//   for (i = 0; i < 10; i++) {
//     for (j = 0; j < 8; j++) {
//       if (collections[i][j] != 0) {
// 	printf("Trying to set up an iterator (%i, %i) \n", i, j);
// 	iter = collections[i][j]->iterator();
// 	if (iter->hasNext()) printf("There is something inside! \n");
// 	while(iter->hasNext()) {
// 	  printf("%i \n", iter->next()->M);
// 	}
//       }
//     }
//   }
//   printf("%i sets, collection[0]->M = %i \n", num_sets, collection[0]->M);
  closedir(root);
  free(str);
}

FeatureCollection::Iterator* StegoModel::getFeatureIterator(int video_birate, int pair, int method, int accept) {
//   if (method < 0 || method >= 10) return 0;
//   if (accept < 0 || accept >= 8)  return 0;
  if (collections[method][accept] != 0) {
    return collections[method][accept]->iterator();
  }
  return 0;
}

featureSet* StegoModel::getCleanSet() {
  return cleanSet;
}


int** StegoModel::getRanges() {
  return ranges;
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

int StegoModel::getSigmaDim() {
  if (mc == 0) return -1;
  return mc->cache;
}

double* StegoModel::getSigma() {
//   if (steg->features == NULL) 
//     return NULL;
//   return steg->features->gauss->qr;
  if (mc == 0) return 0;
  return mc->results;
}

double* StegoModel::getDiag() {
  if (steg->features == NULL) 
    return NULL;
  return steg->features->gauss->qr_diag;
}

// StegoModel::Iterator* StegoModel::getIterator() {
//   return new Iterator(this);
// }
// 
// StegoModel::Iterator::Iterator(StegoModel* model) {
//   level = 0;
//   level1 = -1;
//   level3 = -1;
// //   level0iter = model->collections.iterator();
//   next();
// }
// 
// int StegoModel::Iterator::getLevel0() {
// //   return *level0iter->first;
// }
// 
// int StegoModel::Iterator::getLevel1() {
//   return level1;
// }
// 
// int StegoModel::Iterator::getLevel2() {
// //   return level2iter->first;
// }
// 
// int StegoModel::Iterator::getLevel3() {
//   return level3;
// }
// 
// bool StegoModel::Iterator::hasNext() {
//   return true;
// }
// 
// void StegoModel::Iterator::next() {
//   int i = 3;
//   
// }

// // FeatureCollection::Iterator* StegoModel::Iterator::nextCollectionIterator() {
// //   FeatureCollection::Iterator* iter = level3iter->iterator();
// //   next();
// //   return iter;
// // }
