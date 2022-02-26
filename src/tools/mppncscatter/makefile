# Makefile for the "mppncscatter" program
#
# Written by Hans Vahlenkamp (Hans.Vahlenkamp)
# Geophysical Fluid Dynamics Laboratory / NOAA
# Princeton Forrestal Campus
# Last updated: 05/15/08

#NETCDFPATH = /usr/local/netcdf-4.0.1
SOURCES = mppncscatter.c
OBJECTS = mppncscatter.o
CC = mpicc
CFLAGS =
#LFLAGS = -L${NETCDFPATH}/lib/shared
LIBS = -lnetcdf -lhdf5_hl -lhdf5 -lmpi -lz

mppncscatter: $(OBJECTS)
	$(CC) $(CFLAGS) -o mppncscatter -s $(OBJECTS) $(LFLAGS) $(LIBS)
