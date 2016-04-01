# Makefile for Anaquin, a statistical library for spike-in sequins.
#
#    Copyright (C) 2015- Garvan Institute of Medical Research
#
#    Author: Ted Wong <t.wong@garvan.org.au>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

BOOST = /usr/include/boost

# Required for unit-testing
CATCH = /usr/include/catch/include

# Linear-algebra library
EIGEN = /usr/include/eigen

# Statistics library
SS = ../SS

# Library for SAM/BAM
HLIB = src/htslib

# Include for data model 
HDF5 = /usr/include/hdf5/1.8.14/include

# Library for random generator
KLIB = /usr/include

# Library for data model
HDF5L = /usr/local/Cellar/hdf5/1.8.14/lib

# Where the header are stored
INCLUDE = src

CC = g++
CC_FLAGS = -std=c++11

EXEC         = anaquin
SOURCES      = $(wildcard src/kallisto/src/*.cpp src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/TransQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/LadQuin/*.cpp src/FusQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
SOURCES_TEST = $(wildcard tests/dna/*.cpp tests/parsers/*.cpp tests/TransQuin/*.cpp tests/MetaQuin/*.cpp tests/*.cpp)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/*.c src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -g -L $(HDF5L) -lhdf5 -lz -ldl -o $(EXEC)

%.o: %.c
	gcc -c -I $(HLIB) -I $(INCLUDE) -I $(SS) -I $(EIGEN) -I ${BOOST} -I ${CATCH} -I ${KLIB} $< -o $@

%.o: %.cpp
	$(CC) -g -DK_HACK -DBACKWARD_HAS_BFD -c $(CC_FLAGS) -I $(HDF5) -I $(HLIB) -I $(INCLUDE) -I $(SS) -I $(EIGEN) -I ${BOOST} -I ${CATCH} -I ${KLIB} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
