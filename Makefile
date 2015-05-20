GCC=g++
CPPFLAGS=-std=c++11

CATCH_PATH=/usr/include/catch/include
BOOST_PATH=/usr/local/Cellar/boost/1.57.0/include
EIGEN_PATH=/usr/local/Cellar/eigen/3.2.4/include/eigen3
SS_PATH = ~/Sources/SS
HTLIB_PATH = ~/Sources/HT

all:
	g++ -std=c++11 -I $(BOOST_PATH) -I /usr/include/catch/include/ -I /home/bamboo/boost_1_57_0 -l hts -L $(HTLIB_PATH)  -I $(HTLIB_PATH) -I /home/bamboo/eigen -I src -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I $(SS_PATH)  -I /usr/include -o anaquins src/*.cpp src/rna/*.cpp src/dna/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/rna/t_ralign.cpp tests/t_standard.cpp src/meta/*.cpp


