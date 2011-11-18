CC            = gcc
NVCC          = nvcc
CFLAGS        = -m64
NVCCFLAGS     = -arch=sm_20
LINK          = nvcc
INCPATH       = -I.
LFLAGS        = -m64 -arch=sm_20
LIBS          = $(SUBLIBS)  -lpthread -L/opt/cuda/lib64

####### Output directory

OBJECTS_DIR   = ./

SOURCES       = Stegosaurus.c Stegosaurus_IO.c Stegosaurus_Estimate.cu Stegosaurus_KL.cu
OBJECTS       = Stegosaurus.o Stegosaurus_IO.o Stegosaurus_Estimate.cu.o Stegosaurus_KL.cu.o
TARGET        = Stegosaurus

.SUFFIXES: .o .c .cu

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o "$@" "$<"

all: $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

Stegosaurus.o: Stegosaurus.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus.o Stegosaurus.c

Stegosaurus_IO.o: Stegosaurus_IO.c Stegosaurus_IO.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus_IO.o Stegosaurus_IO.c

Stegosaurus_Estimate.cu.o: Stegosaurus_Estimate.cu Stegosaurus_Estimate.h
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o Stegosaurus_Estimate.cu.o Stegosaurus_Estimate.cu

Stegosaurus_KL.cu.o: Stegosaurus_KL.cu Stegosaurus_KL.h
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o Stegosaurus_KL.cu.o Stegosaurus_KL.cu