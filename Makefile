all:
	g++ -std=c++11 -I /home/bamboo/boost_1_57_0 -l hts -L ../HT  -I ../HT -I /home/bamboo/eigen -I src -I ../../eigen -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I ../QS/include  -I /usr/include -I include/catch/include -I ../../boost -I include/tclap-1.2.1/include/ -o anaquins src/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/t_aligner.cpp tests/t_standard_factory.cpp
