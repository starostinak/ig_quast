INCLUDES=-isystem ext/include/seqan-library-2.0.0/include/ -isystem src/include 
CXXFLAGS=-std=c++11 -Wall -Wextra
CXXFLAGS+=-O3
CXXLIBS=-lpthread -lboost_program_options -lz -lbz2
CXXLIBS+=-lprofiler
CXXLIBS+=-ltcmalloc
CXXLIBS+=-fopenmp
# CXXLIBS+=-fopenmp=libiomp5
# CXXLIBS+=-ljemalloc
# CXX=clang++-3.6
CXX=g++

all: repertoire_evaluator ig_matcher ig_kplus_vj_finder

repertoire_evaluator: src/ig_tools/exact_repertoire_evaluator/main_evaluate_repertoire.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} src/ig_tools/exact_repertoire_evaluator/main_evaluate_repertoire.cpp ${CXXLIBS} -o bin/repertoire_evaluator

ig_matcher: src/ig_tools/ig_matcher/ig_matcher.cpp src/ig_tools/ig_matcher/ig_matcher.hpp src/ig_tools/ig_matcher/banded_half_smith_waterman.hpp
	${CXX} ${CXXFLAGS} ${INCLUDES} src/ig_tools/ig_matcher/ig_matcher.cpp ${CXXLIBS} -o bin/ig_matcher

ig_kplus_vj_finder: src/ig_tools/ig_matcher/ig_kplus_vj_finder.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} src/ig_tools/ig_matcher/ig_kplus_vj_finder.cpp ${CXXLIBS} -o bin/ig_kplus_vj_finder

rep_comp:
	g++ -std=c++11 -g -O0 -DEBUG -Isrc/include -Iext/include/seqan-library-2.0.0/include src/ig_tools/repertoire_comparer/main.cpp src/mph_index/MurmurHash3.cpp -o bin/repertoire_comparer -lboost_thread -lboost_system

clean:
	-rm *.o
	-rm bin/ig_tools/repertoire_evaluator bin/ig_matcher bin/ig_kplus_vj_finder

