# Boost C++ library
BOOST = /usr/local/include/boost_1_64_0

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

HTSLIB = /Users/tedwong/Sources/QA/htslib

# Where the header are
INCLUDE = src

CC = g++
CC_FLAGS = -std=c++11

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -g -lpthread -lz -lhts -L $(HTSLIB) -o $(EXEC)

%.o: %.c
	gcc -c -I src/htslib -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

%.o: %.cpp
	$(CC) -g -DK_HACK -c $(CC_FLAGS) -I src/htslib -I src/stats -I $(INCLUDE) -I $(EIGEN) -I $(HTSLIB) -I ${BOOST} $< -o $@
	#$(CC) -g -DK_HACK -DBACKWARD_HAS_BFD -c $(CC_FLAGS) -I src/htslib -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
