#
#    Copyright (C) 2016- Garvan Institute of Medical Research
#
#    Author: Ted Wong <t.wong@garvan.org.au>
#

#
# To compile the source code, you'll need to have the following:
#
#   1. Boost:  http://www.boost.org/
#   2. Catch:  https://github.com/philsquared/Catch
#   3. Eigen:  https://eigen.tuxfamily.org
#   4. Klib:   https://github.com/attractivechaos/klib
#
# Please email t.wong@garvan.org.au if you have any problems.
#

BOOST = /usr/include/boost

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# Library for random generator
KLIB = /Applications

# Where the header are
INCLUDE = src

CC = g++
CC_FLAGS = -std=c++11

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/*.c src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -g -lpthread -lz -ldl -o $(EXEC)

%.o: %.c
	gcc -c -I src/htslib -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} -I ${KLIB} $< -o $@

%.o: %.cpp
	$(CC) -g -DK_HACK -DBACKWARD_HAS_BFD -c $(CC_FLAGS) -I src/htslib -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} -I ${KLIB} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
