#!/bin/bash

# Copyright 2014 Georgios Karagiannis
#
# Georgios Karagiannis
#
# Postdoctoral research associate
# Department of Mathematics, Purdue University
# 150 N. University Street
# West Lafayette, IN 47907-2067, USA
#
# Telephone: +1 765 494-3405
#
# Email: gkaragia@purdue.edu
#
# Contact email: georgios.stats@gmail.com

CC=icc
CFLAGS=-O2
LDFLAGS=
CPPFLAGS=

FUN=cost_ising2D.c

SOURCES=pisaa.c \
			Crossover_int_operations.c \
			Mutation_int_operations.c \
			MH_int_updates.c \
			Self_adjastment_prosedure.c \
			$(FUN) \
			uniformrng.c \
			nrutil.c
	
OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=exe

# BUILD

build: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# CLEAN

clean:
	rm -rf *o exe
	
# DETAILS

details:
	@echo CC       : $(CC)
	@echo CFLAGS   : $(CFLAGS)
	@echo LDFLAGS  : $(LDFLAGS)
	@echo FUN      : $(FUN)
	@echo CPPFLAGS : $(CPPFLAGS)

# RUN

run:
	./exe

	
	
