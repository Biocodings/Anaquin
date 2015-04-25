GCC=g++
CPPFLAGS=-std=c++11

CATCH_PATH=/usr/include/catch/include
BOOST_PATH=/usr/local/Cellar/boost/1.57.0/include
EIGEN_PATH=/usr/local/Cellar/eigen/3.2.4/include/eigen3

all:
	g++ -std=c++11 -I $(BOOST_PATH) -I /usr/include/catch/include/ -I /home/bamboo/boost_1_57_0 -l hts -L ../HT  -I ../HT -I /home/bamboo/eigen -I src -I ../../eigen -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I ../SS  -I /usr/include -I ../../boost -o anaquins src/*.cpp src/rna/*.cpp src/dna/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/rna/t_raligner.cpp tests/t_standard.cpp src/meta/*.cpp


