all:
	g++ -std=c++11 -I /usr/include -I boost -I libs/tclap-1.2.1/include/ src/*.cpp
