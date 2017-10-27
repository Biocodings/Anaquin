#
# Please modify only BOOST, EIGEN and HTSLIB. You should be able to leave all other options intact. C++ compiler with C++11 support is mandatory.
#

# Boost C++ library
BOOST = /usr/local/include/boost_1_64_0

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# HTSLIB library for reading BAM files
HTSLIB = /Users/tedwong/Sources/QA/htslib

#
# Backward-cpp (https://github.com/bombela/backward-cpp) is useful for C++ stack tracing. Optional.
#

DBACKWARD = #-DBACKWARD_HAS_BFD

CXX      = g++
CC       = $(CXX)
CPPFLAGS = -g -O2 -I src -I src/stats -I src/kallisto
CFLAGS   = -g -O2
CXXFLAGS = -c -std=c++11 
LIBS     = -lpthread -lz -lhts
DFLAGS   = $(DBACKWARD)

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/kallisto/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_LIB)
	$(CXX) $(OBJECTS) $(OBJECTS_LIB) $(CFLAGS) $(DFLAGS) $(LIBS) -L $(HTSLIB) -o $(EXEC)

%.o: %.c
	$(CC)  $(DFLAGS) $(CFLAGS) $(CXXFLAGS) -I $(EIGEN) -I ${BOOST} $< -o $@

%.o: %.cpp
	$(CXX) $(DFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I $(HTSLIB) -I $(EIGEN) -I ${BOOST} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS)