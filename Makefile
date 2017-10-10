# Boost C++ library
BOOST = /usr/local/include/boost_1_64_0

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# HTSLIB library for BAM files
HTSLIB = /Users/tedwong/Sources/QA/htslib

# Where the header are
INCLUDE = src

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	g++ $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) -DBACKWARD_HAS_BFD -g -lpthread -lz -lhts -L $(HTSLIB) -o $(EXEC)

%.o: %.c
	gcc -g -c -DBACKWARD_HAS_BFD -I src/htslib -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

%.o: %.cpp
	g++ -g -DK_HACK -DBACKWARD_HAS_BFD -c -std=c++11 -I src/htslib -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
