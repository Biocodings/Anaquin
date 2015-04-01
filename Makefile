all:
	g++ -std=c++11 -I ~/eigen/ -I /usr/local/Cellar/eigen/3.2.4/include/eigen3  -I ../QS/include  -I /usr/include -I include/catch/include -I boost -I include/tclap-1.2.1/include/ -o anaquins src/*.cpp
