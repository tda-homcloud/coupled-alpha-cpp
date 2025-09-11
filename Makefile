HEADERS := $(wildcard *.hpp)

test_main: test_main.cpp $(HEADERS)
	g++ -g -std=c++17 -Wall -Ieigen -ICGAL-5.6.1/include -o test_main test_main.cpp -lgmp -lmpfr

.PHONY: setup-dependencies test
test: test_main
	./test_main

setup-dependencies:
	curl -L -O https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1.tar.xz
	tar xf CGAL-5.6.1.tar.xz
