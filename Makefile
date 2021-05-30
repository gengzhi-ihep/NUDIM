# ----- Compiler -------------------------------------------------------------
CC = gcc
MPI = /usr/lib64/openmpi/bin/mpicc

# ----- NFFT/FFTW installation directory -------------------------------------
NFFT_DIR = /usr
FFTW_DIR = /usr
GSL_DIR = /usr

# ----- .h header files path -------------------------------------------------
CFLAGS1 = -O3 -I$(NFFT_DIR)/include
CFLAGS2 = -O3 -I$(FFTW_DIR)/include
CFLAGS3 = -O3 -I$(GSL_DIR)/include

# ----- Library directory ----------------------------------------------------
LDFLAGS1 = -L$(NFFT_DIR)/lib64
LDFLAGS2 = -L$(FFTW_DIR)/lib64
LDFLAGS3 = -L$(GSL_DIR)/lib64


# ----- Libs  ----------------------------------------------------------------
LDLIBS1 = -lnfft3
LDLIBS2 = -lfftw3
LDLIBS3 = -lgsl -lgslcblas -lm

# ----------------------------------------------------------------------------
test1: ./src/grid.c ./src/projection_sim.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) ./src/grid.c ./src/projection_sim.c -o ./bin/projection_sim $(LDFLAGS1) $(LDLIBS1) $(LDFLAGS2) $(LDLIBS2) $(LDFLAGS3) $(LDLIBS3)

test2: ./src/grid.c ./src/denoise.c ./src/reconstruction.c
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) ./src/grid.c ./src/denoise.c ./src/reconstruction.c -o ./bin/reconstruction  $(LDFLAGS1) $(LDLIBS1) $(LDFLAGS2) $(LDLIBS2) $(LDFLAGS3) $(LDLIBS3)

test3: ./src/grid.c ./src/denoise.c ./src/reconstruction-parallel.c
	$(MPI) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) ./src/grid.c ./src/denoise.c ./src/reconstruction-parallel.c -o ./bin/reconstruction_parallel  $(LDFLAGS1) $(LDLIBS1) $(LDFLAGS2) $(LDLIBS2) $(LDFLAGS3) $(LDLIBS3)

all:   test1 test2 test3
# ----------------------------------------------------------------------------
clean:
	rm -f ./bin/projection_sim
	rm -f ./bin/reconstruction
	rm -f ./bin/reconstruction_parallel
