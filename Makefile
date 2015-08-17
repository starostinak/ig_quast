.PHONY: dd rd

igtools:
	g++ src/ig_tools/exact_repertoire_evaluator/main_evaluate_repertoire.cpp -std=c++11 -o bin/ig_tools/repertoire_evaluator
	g++ -std=c++11 -Isrc/include -Iext/include/seqan-library-2.0.0/include src/ig_tools/repertoire_comparer/main.cpp src/ssw/ssw.c src/ssw/ssw_cpp.cpp src/mph_index/MurmurHash3.cpp -o bin/ig_tools/repertoire_comparer -lboost_thread -lboost_system

rep_comp:
	g++ -std=c++11 -g -O0 -DEBUG -Isrc/include -Iext/include/seqan-library-2.0.0/include src/ig_tools/repertoire_comparer/main.cpp src/mph_index/MurmurHash3.cpp -o bin/ig_tools/repertoire_comparer -lboost_thread -lboost_system

clean:
	$(MAKE) -C build/debug clean
	$(MAKE) -C build/release clean

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all

#temporary fix
#all: $(MAKE) -C build/release all
