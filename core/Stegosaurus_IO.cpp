#include "Stegosaurus.h"


int readHeader(FILE *file, featureHeader *header) {
  int i;
  
  header->video_bitrate = 3000;
//   printf("pair: %i \n", header->pair);
  fread(&header->pair, sizeof(char), 1, file);
//   printf("pair: %i \n", header->pair);
  fread(&header->slice_type, sizeof(char), 1, file);
//   printf("slice_type: %i \n", header->slice_type);
  fread(&header->method, sizeof(char), 1, file);
//   printf("method: %i \n", header->method);
  if (header->method != 0) {
    fread(&header->using_rate, sizeof(char), 1, file);
    fread(&header->rate, sizeof(double), 1, file);
//     printf("rate: %f \n", header->rate);
    fread(&header->accept, sizeof(char), 1, file);
//     printf("accept: %i \n", header->accept);
  } else {
    header->rate = -1.;
    header->accept = 0;
  }
  fread(&header->qp_offset, sizeof(char), 1, file);
//   printf("offset: %i \n", header->qp_offset);
  fread(&header->qp_range, sizeof(char), 1, file);
//   printf("range: %i \n", header->qp_range);
  
  if (header->pair) {
    fread(&header->ranges[0][0], sizeof(unsigned char), 1, file);
    fread(&header->ranges[0][1], sizeof(unsigned char), 1, file);
    fread(&header->ranges[1][0], sizeof(unsigned char), 1, file);
    fread(&header->ranges[1][1], sizeof(unsigned char), 1, file);
    fread(&header->ranges[2][0], sizeof(unsigned char), 1, file);
    fread(&header->ranges[2][1], sizeof(unsigned char), 1, file);
  } else {
    for (i = 0; i < 16; i++) fread(&header->ranges[0][0]+i, sizeof(unsigned char), 1, file);
    for (i = 0; i < 4; i++)  fread(&header->ranges[1][0]+i, sizeof(unsigned char), 1, file);
    for (i = 0; i < 15; i++) fread(&header->ranges[2][0]+i, sizeof(unsigned char), 1, file);
  }
  
  return 0;
}

int jumpToVector(featureSet *set, long vecnum) {
  int i;
  long seenvs = 0l;
  
//   printf("doing some jump to vec %i... %i \n", vecnum, set->num_files);
  for (i = 0; seenvs <= vecnum && i < set->num_files; i++) {
    seenvs += set->vsPerFile[i];
  }
  seenvs -= set->vsPerFile[--i];
  if (i == set->num_files) {
//     printf("trying to jump out of range! \n");
    return 1;
  }
//   printf("i = %i, seen %i vectors, will jump to pos %i \n", i, seenvs, vecnum-seenvs);
  set->current_file = i;
  fseek(set->files[set->current_file], set->dataOffset + (vecnum-seenvs)*(long)set->dim, SEEK_SET);
//   printf("Done. \n");
//   fsetpos(set->current_file, 
  return 0;
}

void stegoRewind(featureSet *set) {
  int i;

  for (i = 0; i < set->num_files; i++) {
    fseek(set->files[i], set->dataOffset, SEEK_SET);
  }
  set->current_file = 0;
}

// void storeGaussian(char* path, myGaussian* gauss) {
//   
// }
