#include "Stegosaurus.h"


__global__ void compareMax(int dim, double *current_max, double *new_features) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    if (current_max[idx] < new_features[idx])
      current_max[idx] = new_features[idx];
  }
}

__global__ void subtract(int dim, double *vec, double *mu) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim) {
    vec[idx] -= mu[idx];
  }
}

__global__ void initDArray(double *m, int dim, double val) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx < dim)
    m[idx] = val;
}

__global__ void constructQ(double *q, double *diag, double norm, int count) {
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx == 0) {
    if (q[0] > 0) {
      q[0] += norm;
      diag[count] = -norm;
    } else {
      q[0] -= norm;
      diag[count] = norm;
    }
  }
}

/*int computeQPHistogram(stegoContext *steg, double *mu_g, int qp_range, double *result) {
  int i;
  double alpha = 100./((double) steg->features->dim);
  
  for (i = 0; i < qp_range; i++) {
    cublasDasum_v2(steg->gpu_c->handle, steg->features->dim/qp_range, steg->features->mu_g+i*steg->features->dim/qp_range, 1, result+i);
  }
  cublasDscal_v2(steg->gpu_c->handle, 20, &alpha, result, 1);
  
  return 0;
}*/

int estimateMu(stegoContext *steg) {
  int i;
  int read;
//   double sum;
  featureSet *fs = steg->features;
  int dim = fs->dim;
  int tpb = steg->gpu_c->threads_per_block;
//   double zero = 0.;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
  double *current_feature;
  double *mu_g, *max_g, *vec_g;
  
  if (fs->gauss->mu != NULL) return -1;
  
  CUDA_CALL( cudaHostAlloc(&current_feature, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&(fs->gauss->mu), dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&(fs->max_vec), dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&(fs->qp_vec), QP_RANGE*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mu_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&qp_g, 20*sizeof(double)));
  CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&max_g, dim*sizeof(double)));
  
//   fs->gauss->mu = (double*) malloc(dim*sizeof(double));
  for (i = 0; i < dim; i++)
    current_feature[i] = 0.;

  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, mu_g, 1));
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, max_g, 1));

//   printf("uploaded stuff \n");
  for (i = 0; i < fs->M; i++) {
     read = readVector(steg, steg->features, current_feature);
     if (read != dim) printf("read something wrong: %i \n", read);
//      else printf("Read something right! \n");
     CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
     cublasDaxpy(*handle, dim, &(fs->divM), vec_g, 1, mu_g, 1);     // put divM on gpu?
     compareMax<<<BLOCKS(dim,tpb),tpb>>>(dim, max_g, vec_g);
  }
//   computeQPHistogram(steg, fs->mu_g, 20, fs->qp_g);
//   CUBLAS_CALL( cublasGetVector(20, sizeof(double), fs->qp_g, 1, fs->qp_vec, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), mu_g, 1, fs->gauss->mu, 1));
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), max_g, 1, fs->max_vec, 1));
  stegoRewind(fs);
  
  for (i = 0; i < QP_RANGE; i++) {
    fs->qp_vec[i] = 0;
    for (int j = 0; j < dim/QP_RANGE; j++) {
      fs->qp_vec[i] += fs->gauss->mu[i*dim/QP_RANGE+j];
    }
    fs->qp_vec[i] *= 100./((double) dim);
  }
  
  CUDA_CALL( cudaFreeHost(current_feature));
  CUDA_CALL( cudaFree(mu_g));
//   CUDA_CALL( cudaFree(qp_g));
  CUDA_CALL( cudaFree(vec_g));
  CUDA_CALL( cudaFree(max_g));
  
//   printf("Done estimating some mu \n");
  
  return 0;
}

int estimateSigma(stegoContext *steg) {
  int i, j;
  int read;
  int blockwidth, posInSigma;   // size of current block on gpu
  featureSet *fs = steg->features;
  myGaussian *gauss = fs->gauss;
//   double minus_one = -1;
  int dim = fs->dim;
  int tpb = steg->gpu_c->threads_per_block;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
  double *current_feature;
  double *mu_g, *vec_g, *sigma_g;
//   double one = 1.;
//   double *ones_g; // used to sum up sigma using dot product
//   double *ones2_g;
//   double checksum; // all entries should sum to 1, lets verify this
  
  if (!(gauss->mu != NULL && gauss->sigma == NULL)) return -1;
  
//   gauss->sigma = (double*) malloc(dim*dim*sizeof(double));
  CUDA_CALL( cudaHostAlloc(&(gauss->sigma), dim*dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&current_feature, dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMalloc(&mu_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&vec_g, dim*sizeof(double)));
  CUDA_CALL( cudaMalloc(&sigma_g, dim*fs->gpu_matrix_width*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&ones_g, dim*fs->gpu_matrix_width*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&ones2_g, dim*fs->gpu_matrix_width*sizeof(double)));
  
//   initDArray<<<BLOCKS(dim*fs->gpu_matrix_width,tpb),tpb>>>(ones_g, dim*fs->gpu_matrix_width, 1.);
//   initDArray<<<BLOCKS(dim*fs->gpu_matrix_width,tpb),tpb>>>(ones2_g, dim*fs->gpu_matrix_width, 1.);
  
  CUBLAS_CALL( cublasSetVector(dim, sizeof(double), gauss->mu, 1, mu_g, 1));
//   CUBLAS_CALL( cublasDscal(*handle, dim, &(minus_one), mu_g, 1));
  
//   for (i = 0; i < dim*dim; i++) {
//     gauss->sigma[i] = 0;
//   }
  
  for (j = 0; j < (dim+fs->gpu_matrix_width-1)/fs->gpu_matrix_width; j++) { // (dim+fs->gpu_matrix_width-1)/fs->gpu_matrix_width
    posInSigma = j*fs->gpu_matrix_width;
    blockwidth = min(fs->gpu_matrix_width, dim-posInSigma);
//     CUBLAS_CALL( cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->sigma+posInSigma*dim, dim, sigma_g, dim)); // add on cpu??
    initDArray<<<BLOCKS(dim*blockwidth,tpb),tpb>>>(sigma_g, dim*blockwidth, 0.);
    printf("about to calc sigma, blockwidth=%i, dim=%i, pos=%i mw=%i\n", blockwidth, dim, posInSigma, fs->gpu_matrix_width);
    for (i = 0; i < fs->M; i++) { // fs->M
      read = readVector(steg, steg->features, current_feature);
      if (read != dim) printf("read something wrong: %i \n", read);
      CUBLAS_CALL( cublasSetVector(dim, sizeof(double), current_feature, 1, vec_g, 1));
      subtract<<<BLOCKS(dim,tpb),tpb>>>(dim, vec_g, mu_g);
      
//       cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->sigma+j*fs->gpu_matrix_width*dim, dim, sigma_g, dim); // add on cpu??
      
//       CUBLAS_CALL( cublasDdot(*handle, dim, vec_g, 1, ones_g, 1, &checksum));
//       printf("vec checksum: %f, %i, %i \n", checksum, BLOCKS(dim,tpb), tpb);
//       CUBLAS_CALL( cublasDdot(*handle, dim*blockwidth, sigma_g, 1, ones_g, 1, &checksum));
//       printf("sigma checksum: %f \n", checksum);
      CUBLAS_CALL( cublasDger(*handle, dim, blockwidth, &(fs->divM), vec_g, 1, vec_g+j*fs->gpu_matrix_width, 1, sigma_g, dim)); // fs->divM
      if (i%1000 == 0) printf("%i \n", i);
//       cublasGetMatrix(dim, blockwidth, sizeof(double), sigma_g, dim, gauss->sigma+j*fs->gpu_matrix_width*dim, dim); // only do this with multiple blocks
    }
    stegoRewind(fs);
    CUBLAS_CALL( cublasGetMatrix(dim, blockwidth, sizeof(double), sigma_g, dim, gauss->sigma+posInSigma*dim, dim)); // only do this with multiple blocks
  }
  
//   CUDA_CALL( cudaFree(ones2_g));
//   CUDA_CALL( cudaFree(ones_g));
  CUDA_CALL( cudaFree(sigma_g));
  CUDA_CALL( cudaFree(mu_g));
  CUDA_CALL( cudaFree(vec_g));
  CUDA_CALL( cudaFreeHost(current_feature));
  
  return 0;
}

klContext* initKLContext(stegoContext *steg) {
  int i;
  klContext *klc = (klContext*) malloc(sizeof(klContext));
  
  klc->dim = steg->features->dim;
  klc->nstr = steg->gpu_c->num_streams;
  klc->streams = (cudaStream_t*) malloc(klc->nstr*sizeof(cudaStream_t));
  klc->tmp_g = (double**) malloc(klc->nstr*sizeof(double*));
  klc->vec_g = (double**) malloc(klc->nstr*sizeof(double*));
  klc->dotp_g = (double**) malloc(klc->nstr*sizeof(double*));
  
  klc->handle = &steg->gpu_c->handle;
  
  for (i = 0; i < klc->nstr; i++) {
    CUDA_CALL( cudaStreamCreate(&(klc->streams[i])));
    CUDA_CALL( cudaMalloc(&(klc->tmp_g[i]), klc->dim*sizeof(double)));
    CUDA_CALL( cudaMalloc(&(klc->vec_g[i]), klc->dim*sizeof(double)));
    CUDA_CALL( cudaMalloc(&(klc->dotp_g[i]), sizeof(double)));
  }
  
  return klc;
}

void closeKLContext(klContext *klc) {
}

/*void constructQ(double *q, double *diag, double norm, int count) {
  if (q[0] > 0) {
    q[0] += norm;
    diag[count] = -norm;
  } else {
    q[0] -= norm;
    diag[count] = norm;
  }
}*/

// apply q to vec and store result in vec. vec doesn't change if q=0 (which actually shouldn't happen!)
void applyQ(cublasHandle_t handle, int dim, double *q, double *vec, double *tmp, double *dotp, double *mintwo) {
  CUBLAS_CALL( cublasDcopy(handle, dim, q, 1, tmp, 1));
  CUBLAS_CALL( cublasDdot(handle, dim, tmp, 1, vec, 1, dotp));
  CUBLAS_CALL( cublasDscal(handle, dim, dotp, tmp, 1)); // /vTv
  CUBLAS_CALL( cublasDaxpy(handle, dim, mintwo, tmp, 1, vec, 1));
}

void applyQs(klContext *klc, double *qs_g, double *vs_g, int numqs, int numvs, int cdim) {
  int j, k, s;
  int current_dim;// = cdim;
  double mintwo = -2.;
  double dotp;
  double *currentq_g;
  double *current_vec;
  
  for (j = 0; j < numqs; j++) {
    printf("doing q num %i \n", j);
    current_dim = cdim - j;
    currentq_g = qs_g + j*(klc->dim+1);
//     CUBLAS_CALL( cublasSetVector(current_dim, sizeof(double), current_mat, 1, currentq_g, 1));
    for (k = 0; k < numvs; k+=klc->nstr) {
      for (s = 0; s < klc->nstr; s++) {
        CUBLAS_CALL( cublasSetStream(*klc->handle, klc->streams[s]));
// 	printf("applyQs: set strem nr. %i \n", s);
        if (k + s < numvs) {
          current_vec = vs_g + (k+s)*klc->dim+j;
//     CUBLAS_CALL( cublasSetVectorAsync(current_dim, sizeof(double), current_vec, 1, vec_g[s], 1, streams[s]));
          applyQ(*klc->handle, current_dim, currentq_g, current_vec, klc->tmp_g[s], &dotp, &mintwo); // klc->dotp_g[s]
        }
      }
    }
  }
  for (j = 0; j < klc->nstr; j++) {
    CUDA_CALL( cudaStreamSynchronize(klc->streams[j]));
  }
  CUBLAS_CALL( cublasSetStream(*klc->handle, NULL));
}

void qrHouseholder(stegoContext *steg) {
  int i, j;// k, s;
  int blockwidth, posInSigma;
//   int nstr = steg->gpu_c->num_streams;
  klContext *klc;
  featureSet *fs = steg->features;
  myGaussian *gauss = fs->gauss;
  int dim = fs->dim;
  int tpb = steg->gpu_c->threads_per_block;
  cublasHandle_t *handle = &(steg->gpu_c->handle);
//   double *current_mat;  // cpu
//   double *current_vec;  // cpu
  double *currentq_g;   // gpu
  double *diag_g;
  double *qr_g;
  double *vectorBlock;
  int     current_dim;
//   double mintwo = -2.;
//   cudaStream_t streams[nstr];
//   double *tmp_g[nstr];
//   double *vec_g[nstr];
//   double *dotp_g[nstr];
  double result;
  
  if (gauss->sigma == NULL) estimateSigma(steg);
  
  // init
  klc = initKLContext(steg);
  CUDA_CALL( cudaHostAlloc(&(gauss->qr), dim*dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaHostAlloc(&(gauss->qr_diag), dim*sizeof(double), cudaHostAllocDefault));
  CUDA_CALL( cudaMemcpy(gauss->qr, gauss->sigma, dim*dim*sizeof(double), cudaMemcpyHostToHost)); // just overwrite sigma??
  CUDA_CALL( cudaMalloc(&qr_g, dim*steg->features->gpu_matrix_width*sizeof(double)));
  CUDA_CALL( cudaMalloc(&vectorBlock, dim*steg->features->gpu_matrix_width*sizeof(double)));
  CUDA_CALL( cudaMalloc(&diag_g, dim*sizeof(double)));
//   CUDA_CALL( cudaMalloc(&currentq_g, dim*sizeof(double)));
  
  initDArray<<<BLOCKS(dim,tpb),tpb>>>(diag_g, dim, 0.);
  
  printf("qr: first for loop \n");
  for (i = 0; i < (dim+fs->gpu_matrix_width-1)/fs->gpu_matrix_width; i++) {
    posInSigma = i*fs->gpu_matrix_width;
    blockwidth = min(fs->gpu_matrix_width, dim-posInSigma);
    CUBLAS_CALL( cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->qr+posInSigma*dim, dim, qr_g, dim)); 
    // apply previous q's
    printf("qr: applying prev q's \n");
    for (j = 0; j < i; j++) {
      printf("j = %i ", j);
      CUBLAS_CALL( cublasSetMatrix(dim, fs->gpu_matrix_width, sizeof(double), gauss->qr+j*fs->gpu_matrix_width*dim, dim, vectorBlock, dim)); // only do this with multiple blocks
      printf("applying q's (%i) \n", fs->gpu_matrix_width);
      applyQs(klc, vectorBlock+j*fs->gpu_matrix_width, qr_g+j*fs->gpu_matrix_width, fs->gpu_matrix_width, blockwidth, dim-j*fs->gpu_matrix_width);
    }
    
    //// DEBUG
//     CUBLAS_CALL( cublasSetMatrix(dim, blockwidth, sizeof(double), gauss->qr+posInSigma*dim, dim, vectorBlock, dim));
    ////
    
    printf("qr: start proper work round %i \n", i);
    for (j = 0; j < blockwidth; j++) { // blockwidth
      current_dim = dim - posInSigma - j;
      currentq_g = qr_g + (j+1)*dim - current_dim;
      // find next q
      printf("qr: current_dim = %i, posInSigma = %i, j= %i \n", current_dim, posInSigma, j);
      CUBLAS_CALL( cublasDnrm2(*handle, current_dim, currentq_g, 1, &result)); // klc->dotp_g[0]
      printf("norm = %g \n", result);
//       printPointerInfo(klc->dotp_g[0]);
      if (result > 0.00000000000001) { // != 0.
	printf("norm not 0 \n");
        constructQ<<<1,1>>>(currentq_g, diag_g, result, posInSigma+j); // klc->dotp_g[0]
	CUBLAS_CALL( cublasDnrm2(*handle, current_dim, currentq_g, 1, &result)); // klc->dotp_g[0]
	result = 1./result;
	CUBLAS_CALL( cublasDscal(*handle, current_dim, &result, currentq_g, 1));
        applyQs(klc, currentq_g, currentq_g+dim, 1, blockwidth-j-1, current_dim);
      } else {
	initDArray<<<1,1>>>(diag_g+posInSigma+j, 1, 0.);
      }
      //// DEBUG
//       applyQs(klc, currentq_g, vectorBlock+dim-current_dim, 1, blockwidth, current_dim);
      ////
    }
    printf("done proper work. \n");
    CUBLAS_CALL( cublasGetMatrix(dim, blockwidth, sizeof(double), qr_g, dim, gauss->qr+posInSigma*dim, dim)); 
//     CUBLAS_CALL( cublasGetMatrix(dim, blockwidth, sizeof(double), vectorBlock, dim, gauss->qr+posInSigma*dim, dim)); 
  }
  printf("done with first loop \n");
  CUBLAS_CALL( cublasGetVector(dim, sizeof(double), diag_g, 1, gauss->qr_diag, 1)); // diag
//   printf("First diag element: %f \n", gauss->qr_diag[0]);
}

void qrUnitTest() {
  printf("doing unit test \n");
  int i, j;
  double mu[] = {1., 1., 1., 1., 1., 1.};
  double sigma[] = {1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1,  1, 2, 3, 4, 5, 6,  1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1,  1, 2, 3, 4, 5, 6};
  stegoContext *steg = init_stego();
  myGaussian gauss;
  featureSet set;
  
  gauss.mu = mu;
  gauss.sigma = sigma;
  set.dim = 6;
  set.gauss = &gauss;
  steg->features = &set;
  steg->features->gpu_matrix_width = 6;
  
  qrHouseholder(steg);
  
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      printf("%f ", gauss.qr[6*j + i]);
    }
    printf("\n");
  }
}
