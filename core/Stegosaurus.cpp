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

int FeatureCollection::addFeatureFile(const char *path, featureHeader *header) {
  int bin = (int) (header->rate/BIN_WIDTH + 0.5);
  printf("rate = %g, bin = %i \n", header->rate, bin);
  featureSet *set;

  printf("adding feature file \n");
  if (collection[bin] == 0) {
    printf("no collection for this yet! \n");
    set = openFeatureSet(path);
    set->rate = header->rate;
    set->id = bin;
    collection[bin] = set;
    printf("just defined collection[%i] \n", bin);
  } else {
    printf("There already is something \n");
//     newFeatureFile(collection[bin], path);
  }
  printf("added new feature file \n");
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
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 8; j++) {
      collections[i][j] = 0;
    }
  }
}

StegoModel::~StegoModel() {
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
      if (collections[header.method][header.accept] == 0) {
	collections[header.method][header.accept] = new FeatureCollection(&header);
	printf("Created new collection for method %i and accept %i \n", header.method, header.accept);
      }
      collections[header.method][header.accept]->addFeatureFile(str, &header);
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

FeatureCollection::Iterator* StegoModel::getFeatureIterator(int method, int accept) {
  if (method < 0 || method >= 10) return 0;
  if (accept < 0 || accept >= 8)  return 0;
  if (collections[method][accept] != 0) {
    return collections[method][accept]->iterator();
  }
  return 0;
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

