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
}

FeatureCollection::~FeatureCollection() {
  int i;
  map< int, featureSet* >::iterator fiter;
  
  for (fiter = collection.begin(); fiter != collection.end(); fiter++) {
    closeFeatureSet(fiter->second);
  }
}

int FeatureCollection::addFeatureFile(const char* path, featureHeader* header, stegoContext* steg, featureSet* cleanSet) {
  int bin;
  featureSet *set;

  bin = (int) (header->prob*10000 + 0.5);
//   printf("using bin %i \n", bin);
  if (collection[bin] == 0) {
    if (cleanSet->header->slice_type == header->slice_type) { // maybe want more chekcs here
      set = openFeatureSet(path, steg);
      set->mask_vec = cleanSet->mask_vec;
      set->mask_vec = cleanSet->mask_vec;
      set->mask_counts = cleanSet->mask_counts;
      set->max_g    = cleanSet->max_g;
      set->min_g    = cleanSet->min_g;
      set->mu_g     = cleanSet->mu_g;
      set->mu_vec   = cleanSet->mu_vec;
      set->var_g    = cleanSet->var_g;
      
      set->prob = header->prob;
      set->id = bin;
      collection[bin] = set;
    }
  } else {
    if (collection[bin]->header->slice_type == header->slice_type)
      newFeatureFile(steg, collection[bin], path);
  }
}

int FeatureCollection::getNumSets() {
   return collection.size();
}

featureSet* FeatureCollection::getFeatureSet(int index) {
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
  seenPaths = new set<string>();
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
  map< pair< int, int >, FeatureCollection* >::iterator fiter;
  FeatureCollection::Iterator *citer;
//   featureSet *stego;
  steg->features = cleanSet;
  startAction(steg->features);
  estimateMu(steg);
  estimateSigma(steg);
//   qrHouseholder(steg);
  endAction(steg->features);
  for (fiter = collections.begin(); fiter != collections.end(); fiter++) {
    if (fiter->second != 0) {
//       printf("<%i, %i> \n", fiter->first.first, fiter->first.second);
      citer = fiter->second->iterator();
      while (citer->hasNext()) {
	steg->features = citer->next();
// 	printf("About to estimate mu with dim=%i, M=%i \n", steg->features->dim, steg->features->M);
	startAction(steg->features);
	estimateMu(steg);
	endAction(steg->features);
// 	progressChanged((double) i / (double) j);
      }
    }
  }
  steg->features = cleanSet;
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

void StegoModel::doMMDs() {
  map< pair< int, int >, FeatureCollection* >::iterator fiter;
  FeatureCollection::Iterator *citer;
//   featureSet *stego;
  
  mc = (mmdContext*) malloc(sizeof(mmdContext));
  mc->clean = cleanSet;
  startAction(mc->clean);
  initMMD(steg, *mc);
  estimateGamma(steg, *mc);
    
  for (fiter = collections.begin(); fiter != collections.end(); fiter++) {
    if (fiter->second != 0) {
      printf("<%i, %i> \n", fiter->first.first, fiter->first.second);
      citer = fiter->second->iterator();
      while (citer->hasNext()) {
	mc->stego = citer->next();
	printf("doing set %g \n", mc->stego->header->prob);
	startAction(mc->stego);
	estimateMMD(steg, *mc);
	endAction(mc->stego);
// 	doMMD(cleanSet, stego); // don't want to call this here, should be more low-level
      }
    }
  }
  endAction(mc->clean);
  closeMMD(*mc);
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
//     printf("updateing some view \n");
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
      openFile(str, i, num_sets);
//       file = fopen(str, "r");
//       readHeader(file, &header);
//       if (ranges == 0) {
// 	printf("first file is being added! \n");
// 	ranges = (int**) malloc(3*sizeof(int*));
// 	for (k = 0; k < 3; k++) {
// 	  ranges[k] = (int*) malloc(num_coefs[k]*sizeof(int));
// 	  for (j = 0; j < num_coefs[k]; j++)
// 	    ranges[k][j] = header.ranges[k][j];
// 	}
//       }
//       fclose(file);
// 
// //      printf("method: %i \n", header.method);
// //       printf("qp range: %i \n", header.qp_range);
//       if (header.method == 0) {
// 	if (cleanSet == 0) {
// 	  cleanSet = openFeatureSet(str, steg);
// 	} else {
// // 	  printf("adding new feature file to clean set \n");
// 	  newFeatureFile(steg, cleanSet, str);
// 	}
//       } else {
// 	if (collections[pair<int, int>(header.method,header.accept)] == 0) {
// 	  collections[pair<int, int>(header.method,header.accept)] = new FeatureCollection(&header);
// 	  printf("Created new collection for method %i and accept %i \n", header.method, header.accept);
// 	}
//         collections[pair<int, int>(header.method, header.accept)]->addFeatureFile(str, &header, steg, cleanSet);
//       }
//       i++;
//       progressChanged((double) i / (double) num_sets);
    }
  }
  collectionChanged();

  closedir(root);
  free(str);
}

void StegoModel::openFile(const char* path, int i, int num_sets) {
  int j, k;
  FILE *file;
  featureHeader header;
  
  if (seenPaths->find(string(path)) != seenPaths->end()) {
//     printf("Das kenne ich doch schon %s :-@ \n", path);
    return;
  }
  
  file = fopen(path, "r");
  readHeader(file, &header);
  if (ranges == 0) {
    ranges = (int**) malloc(3*sizeof(int*));
    for (k = 0; k < 3; k++) {
      ranges[k] = (int*) malloc(num_coefs[k]*sizeof(int));
      for (j = 0; j < num_coefs[k]; j++)
        ranges[k][j] = header.ranges[k][j];
    }
    for (k = 0; k < 3; k++) {
      for (j = 0; j < num_coefs[k]; j++) {
	ranges[k][j] = 2 * (int) header.ranges[k][j] + 1;
// 	printf("ranges[%i][%i] = %i \n", k, j, ranges[k][j]);
      }
    }
  }
  fclose(file);

  if (header.method == 0) {
    if (cleanSet == 0) {
      cleanSet = openFeatureSet(path, steg);
    } else {
      newFeatureFile(steg, cleanSet, path);
    }
  } else {
    if (collections[pair<int, int>(header.method,header.accept)] == 0) {
      collections[pair<int, int>(header.method,header.accept)] = new FeatureCollection(&header);
    }
    collections[pair<int, int>(header.method, header.accept)]->addFeatureFile(path, &header, steg, cleanSet);
  }
  seenPaths->insert(string(path));
  progressChanged((double) i / (double) num_sets);
  collectionChanged(); // maybe a bit inefficient to run this on every file, might be hundreds
}

StegoModel::Iterator::Iterator(StegoModel* m) {
  model = m;
  iter = m->collections.begin();
}

bool StegoModel::Iterator::hasNext() {
  return (iter != model->collections.end());
}

FeatureCollection* StegoModel::Iterator::next() {
  FeatureCollection *fc = iter->second;
  x = iter->first;
  iter++;
  return fc;
}

std::pair< int, int > StegoModel::Iterator::getX() {
  return x;
}

StegoModel::Iterator* StegoModel::iterator() {
  return new StegoModel::Iterator(this);
}

// FeatureCollection::Iterator* StegoModel::getFeatureIterator(int video_birate, int ppair, int method, int accept) {
// //   if (method < 0 || method >= 10) return 0;
// //   if (accept < 0 || accept >= 8)  return 0;
// //   if (collections[method][accept] != 0) {
// //     return collections[method][accept]->iterator();
// //   }
//   if (collections[pair<int, int>(method, accept)] != 0) {
//     return collections[pair<int, int>(method, accept)]->iterator();
//   }
//   return 0;
// }

featureSet* StegoModel::getCleanSet() {
  return cleanSet;
}

int StegoModel::getQPRange() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->header->qp_range;
}


int** StegoModel::getRanges() {
  return ranges;
}


int StegoModel::getDimension() {
  if (steg->features == NULL)
    return -1;
  return steg->features->dim;
}

int StegoModel::getHistDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->hist_dim;
}

int StegoModel::getPairDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->pair_dim;
}

int StegoModel::getUvsVDim() {
  if (cleanSet == 0)
    return -1;
  return cleanSet->uvsv_dim;
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

// FeatureCollection* StegoModel::getCollection() {
// //   return fcol; // fcol
//   return 0;
// }

int StegoModel::getSigmaDim() {
  if (mc == 0) return -1;
  return mc->cache;
}

double* StegoModel::getSigma() {
  if (steg->features == 0) 
    return 0;
  return steg->features->gauss->qr;
//   if (mc == 0) return 0;
//   return mc->results;
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
