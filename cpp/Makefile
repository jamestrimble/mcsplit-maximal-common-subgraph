CXX := g++
CXXFLAGS := -O3
#CXXFLAGS := -g -O0 -fsanitize=address -static-libasan
all: mcsp mcsp_implicit

mcsp: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp graph.c mcsp.c -pthread

mcsp_implicit: mcsp_implicit.c graph_implicit.c graph_implicit.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_implicit graph_implicit.c mcsp_implicit.c -pthread
