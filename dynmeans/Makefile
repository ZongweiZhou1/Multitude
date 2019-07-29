_dynmeans.so : dynmeans.o dynmeans_wrap.o
	gcc -shared dynmeans.o dynmeans_wrap.o -o _dynmeans.so -lstdc++ -lpython3.5m

dynmeans.o : dynmeans.cpp
	g++ -c -fPIC dynmeans.cpp -I/usr/include/python3.5

dynmeans_wrap.o : dynmeans_wrap.cxx
	g++ -c -fPIC dynmeans_wrap.cxx -I/usr/include/python3.5
dynmeans_wrap.cxx dynmeans.py : dynmeans.i dynmeans.hpp
	swig -c++ -python dynmeans.i

clean:
	rm -f *.o *.so dynmeans_wrap.* dynmeans.py*

test:
	python test.py

.DEFAULT_GOAL := _dynmeans.so

.PHONY: clean test