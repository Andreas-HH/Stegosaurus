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

int newFeatureFile(featureSet* set, const char* path) {
  int i;
  int dim = 0;
  int read;
  int vec[set->dim];
  long localM = 0l;
  FILE *file;// = fopen(path, "r");
  featureHeader header;
    
  if (set->num_files == MAX_FILES) return -1;
  file = fopen(path,"r");
  readHeader(file, &header);
//   // double-check dimension
//   if (header.pair) {
//     dim = (2*header.ranges[0][0]+1)*(2*header.ranges[0][1]+1) + 
//           (2*header.ranges[1][0]+1)*(2*header.ranges[1][1]+1) + 
// 	  (2*header.ranges[2][0]+1)*(2*header.ranges[2][1]+1);
//   } else {
//     for (i = 0; i < 16; i++) dim += 2*header.ranges[0][i]+1;
//     for (i = 0; i < 4; i++)  dim += 2*header.ranges[0][i]+1;
//     for (i = 0; i < 15; i++) dim += 2*header.ranges[0][i]+1;
//   }
//   dim *= header.qp_range;
//   if (dim != set->dim) {
//     printf("Dimension mismatch! \n");
//     fclose(file);
//     return -2;
//   }
  
  while ((read = fread(vec, sizeof(int), dim, file))>0) // && M<10000
    localM++;
  set->M += (int) localM;
  set->vsPerFile[set->num_files] = localM;
  set->divM = 1./set->M;
  rewind(file);
  readHeader(file, &header);
  
  set->files[set->num_files] = file;
  set->num_files++;
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
//   featureHeader header;
//   
  for (i = 0; i < set->num_files; i++) {
//     fsetpos(set->files[i], &set->dataOffset);
    fseek(set->files[i], set->dataOffset, SEEK_SET);
//     rewind(set->files[i]);
//     readHeader(set->files[i], &header);
  }
//   set->current_file = 0;
}

// void storeGaussian(char* path, myGaussian* gauss) {
//   
// }
