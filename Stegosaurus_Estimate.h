#ifndef STEGOSAURUS_ESTIMATE
#define STEGOSAURUS_ESTIMATE

#include<Stegosaurus_IO.h>

int estimateGaussian(char *features);

extern "C" {  
  int init_gpu();
}

#endif /* STEGOSAURUS_ESTIMATE */