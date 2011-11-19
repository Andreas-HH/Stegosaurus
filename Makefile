CC            = gcc
NVCC          = nvcc
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT
NVCCFLAGS     = -arch=sm_20
LINK          = nvcc
INCPATH       = -I.
LFLAGS        = -arch=sm_20
LIBS          = -lpthread -L/opt/cuda/lib64 -lcublas

####### Output directory

OBJECTS_DIR   = ./

SOURCES       = Stegosaurus.c Stegosaurus_IO.c Stegosaurus_Estimate.cu Stegosaurus_KL.cu
OBJECTS       = Stegosaurus.o Stegosaurus_Estimate.cu.o Stegosaurus_Estimate.o Stegosaurus_IO.o Stegosaurus_KL.cu.o Stegosaurus_KL.o
TARGET        = Stegosaurus

.SUFFIXES: .o .c .cu

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o "$@" "$<"

all: $(TARGET)

$(TARGET): $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

Stegosaurus.o: Stegosaurus.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus.o Stegosaurus.c

Stegosaurus_IO.o: Stegosaurus_IO.c Stegosaurus_IO.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus_IO.o Stegosaurus_IO.c

Stegosaurus_Estimate.o: Stegosaurus_Estimate.c Stegosaurus_Estimate.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus_Estimate.o Stegosaurus_Estimate.c

Stegosaurus_KL.o: Stegosaurus_KL.c Stegosaurus_KL.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o Stegosaurus_KL.o Stegosaurus_KL.c

Stegosaurus_Estimate.cu.o: Stegosaurus_Estimate.cu Stegosaurus_Estimate.h
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o Stegosaurus_Estimate.cu.o Stegosaurus_Estimate.cu

Stegosaurus_KL.cu.o: Stegosaurus_KL.cu Stegosaurus_KL.h
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o Stegosaurus_KL.cu.o Stegosaurus_KL.cu