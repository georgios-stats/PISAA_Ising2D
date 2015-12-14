#!/bin/bash

# --------------------------------------------------------------------------------
# 
# Copyrigtht 2014 Georgios Karagiannis
# 
# This file is part of PISAA_Ising2D.
# 
# PISAA_Ising2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# PISAA_Ising2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PISAA_Ising2D.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------

# Georgios Karagiannis
#
# Postdoctoral research associate
# Department of Mathematics, Purdue University
# 150 N. University Street
# West Lafayette, IN 47907-2067, USA
#
# Telephone: +1 (765) 496-1007
#
# Email: gkaragia@purdue.edu
#
# Contact email: georgios.stats@gmail.com
#

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

	
	
