#
#    Copyright (C) 2016- Garvan Institute of Medical Research
#
#    Author: Ted Wong <t.wong@garvan.org.au>
#

#
# https://s3.amazonaws.com/sequins/software/CompileAnaquin.pdf has the instructions
#

BOOST = /usr/local/Cellar/boost/1.60.0_1/include

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# Library for VCF parsing
VCFLIB = /Users/tedwong/Sources/VCF/vcflib

# Where the header are
INCLUDE = src

CC = g++
CC_FLAGS = -std=c++11

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/*.c src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -L $(VCFLIB) -g -lpthread -lz -lvcflib -o $(EXEC)

%.o: %.c
	gcc -c -I src/htslib -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

%.o: %.cpp
	$(CC) -g -DK_HACK -c $(CC_FLAGS) -I src/htslib -I $(VCFLIB)/include -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@
	#$(CC) -g -DK_HACK -DBACKWARD_HAS_BFD -c $(CC_FLAGS) -I src/htslib -I $(VCFLIB)/include -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@
	
clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
