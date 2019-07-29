_maxmatching.so : maxmatching.o maxmatching_wrap.o
	gcc -shared maxmatching.o maxmatching_wrap.o -o _maxmatching.so -lstdc++ -lpython3.5m -llpsolve55

maxmatching.o : maxmatching.cpp
	g++ -c -fPIC maxmatching.cpp -I/usr/include/python3.5

maxmatching_wrap.o : maxmatching_wrap.cxx
	g++ -c -fPIC maxmatching_wrap.cxx -I/usr/include/python3.5

maxmatching_wrap.cxx maxmatching.py : maxmatching.i maxmatching.hpp
	swig -c++ -python maxmatching.i

clean:
	rm -f *.o *.so maxmatching_wrap.* maxmatching.py*

test:
	python test.py

.DEFAULT_GOAL := _maxmatching.so

.PHONY: clean test