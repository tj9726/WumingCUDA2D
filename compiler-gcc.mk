# -*- Makefile -*-

# HDF5
HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR      = ar
CC      = mpicc -fopenmp
FC      = mpif90 -fopenmp
FCFLAGS = -cpp -I$(BASEDIR)/include -I$(HDF5DIR)/include
LDFLAGS = -L$(BASEDIR)/lib -L$(HDF5DIR)/lib
