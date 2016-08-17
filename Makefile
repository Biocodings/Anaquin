#
#    Copyright (C) 2016- Garvan Institute of Medical Research
#
#    Author: Ted Wong <t.wong@garvan.org.au>
#

BOOST = /usr/include/boost

# Required for unit-testing
CATCH = /Applications/catch/include

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# Statistics library
SS = ../SS

# Library for random generator
KLIB = /Applications

# Where the header are stored
INCLUDE = src

CC = g++
CC_FLAGS = -std=c++11

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
SOURCES_TEST = $(wildcard tests/dna/*.cpp tests/parsers/*.cpp tests/RnaQuin/*.cpp tests/*.cpp)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/*.c src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -g -L $(SS) -lnmaths -lz -ldl -o $(EXEC)

%.o: %.c
	gcc -c -I src/htslib -I $(INCLUDE) -I $(SS) -I $(EIGEN) -I ${BOOST} -I ${CATCH} -I ${KLIB} $< -o $@

%.o: %.cpp
	$(CC) -g -DK_HACK -DBACKWARD_HAS_BFD -c $(CC_FLAGS) -I src/htslib -I tests -I $(INCLUDE) -I $(SS) -I $(EIGEN) -I ${BOOST} -I ${CATCH} -I ${KLIB} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
