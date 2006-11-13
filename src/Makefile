# Makefile for pcamal_proto
EXEC= pcamal_proto

INCDIR= -I/usr/include/openmpi
LIBDIR= -L/usr/lib/openmpi -lmpi

SRC= pcamal_proto.c
OBJ= pcamal_proto.o

CFLAGS= -Wall

CC=gcc

.SILENT:

$(EXEC): $(OBJ)
	@echo linking $@ ... 
	$(CC) -o $@ $(OBJ) $(LIBDIR)

%.o: %.c
	@echo compiling $@ ...  
	$(CC) -c $(CFLAGS) $(INCDIR) $< -o $@

clean::
	@echo cleaning ...
	rm -f *.o 
	rm -f $(EXEC)
