GCC=g++
CPPFLAGS=-std=c++11

CATCH_PATH = {CATCH_PATH}
BOOST_PATH = {BOOST_PATH}
EIGEN_PATH = {EIGEN_PATH}
SS_PATH    = {SS_PATH}
HTLIB_PATH = {HTLIB_PATH}

all:
	g++ -std=c++11 -I $(BOOST_PATH) -I /usr/include/catch/include/ -I /home/bamboo/boost_1_57_0 -l hts -L $(HTLIB_PATH)  -I $(HTLIB_PATH) -I /home/bamboo/eigen -I src -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I $(SS_PATH)  -I /usr/include -o anaquins src/*.cpp src/rna/*.cpp src/dna/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/rna/t_ralign.cpp tests/t_standard.cpp src/meta/*.cpp


