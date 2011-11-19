#ifndef STEGOSAURUS_IO
#define STEGOSAURUS_IO

#include<Stegosaurus.h>


void storeGaussian(char *path, MyGaussian *gauss);
FeatureSet* openFeatureSet(char *path);
int closeFeatureSet(FeatureSet *set);
int readVector(FeatureSet *set);

#endif /* STEGOSAURUS_IO */