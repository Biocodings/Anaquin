#
# Please modify only BOOST, EIGEN and HTSLIB. You should be able to leave all other options intact. C++ compiler with C++11 support is mandatory.
#

# Boost C++ library
BOOST = /usr/local/include/boost_1_64_0

# Linear-algebra library
EIGEN = /usr/local/Cellar/eigen/3.2.8/include/eigen3

# HTSLIB library for reading BAM files
HTSLIB = /Users/tedwong/Sources/QA/htslib

CC     = g++
CFLAGS = -g -O2
CPPFLAGS = -c -std=c++11
DFLAGS = 
#DFLAGS = -DBACKWARD_HAS_BFD # https://github.com/bombela/backward-cpp
LIBS   = -lpthread -lz -lhts

# Where the header are (no need to modify this)
INCLUDE = src

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/RnaQuin/*.cpp src/VarQuin/*.cpp src/MetaQuin/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp src/cufflinks/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
SOURCES_LIB  = $(wildcard src/htslib/cram/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CC) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB) $(CFLAGS) $(DFLAGS) $(LIBS) -L $(HTSLIB) -o $(EXEC)

%.o: %.c
	$(CC) $(CFLAGS) -c $(DFLAGS) -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(DFLAGS) $(CPPFLAGS) -I $(HTSLIB) -I src/stats -I $(INCLUDE) -I $(EIGEN) -I ${BOOST} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
