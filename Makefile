all:
	g++ -std=c++11 -I /usr/include -I include/catch/include -I boost -I include/tclap-1.2.1/include/ src/*.cpp
