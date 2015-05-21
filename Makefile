# Makefile for anaquin, a statistical library for spike-in sequins.
#
#    Copyright (C) 2013-2015 Garvan Institute of Medical Research
#
#    Author: Ted Wong <t.wong@garvan.org.au>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

GCC=g++
CPPFLAGS=-std=c++11

CATCH_PATH = {CATCH_PATH}
BOOST_PATH = {BOOST_PATH}
EIGEN_PATH = {EIGEN_PATH}
SS_PATH    = {SS_PATH}
HTLIB_PATH = {HTLIB_PATH}

all:
	g++ -std=c++11 -I $(BOOST_PATH) -I {CATCH_PATH} -I /home/bamboo/boost_1_57_0 -l hts -L $(HTLIB_PATH)  -I $(HTLIB_PATH) -I /home/bamboo/eigen -I src -I ${EIGEN_PATH}  -I $(SS_PATH)  -I /usr/include -o anaquins src/*.cpp src/rna/*.cpp src/dna/*.cpp src/parsers/*.cpp src/writers/*.cpp tests/rna/t_ralign.cpp tests/t_standard.cpp src/meta/*.cpp
