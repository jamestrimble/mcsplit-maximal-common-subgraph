CXX := g++
CXXFLAGS := -O3
#CXXFLAGS := -g -O0 -fsanitize=address -static-libasan
all: mcsp

mcsp: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp graph.c mcsp.c -pthread
