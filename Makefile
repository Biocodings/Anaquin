all:
	g++ -std=c++11 -I /usr/include/catch/include/ -I /home/bamboo/boost_1_57_0 -l hts -L ../HT  -I ../HT -I /home/bamboo/eigen -I src -I ../../eigen -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I ../SS  -I /usr/include -I ../../boost -o anaquins src/*.cpp src/rna/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/rna/t_aligner.cpp tests/t_standard_factory.cpp
