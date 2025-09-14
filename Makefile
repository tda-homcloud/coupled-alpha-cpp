HEADERS := $(wildcard *.hpp)
FLAGS := -g -std=c++17 -Wall -Ieigen -ICGAL-5.6.1/include
test_main: test_main.cpp $(HEADERS)
	g++ $(FLAGS) -o test_main test_main.cpp -lgmp -lmpfr

coupled_alpha: main.cpp $(HEADERS)
	g++ $(FLAGS) -o coupled_alpha main.cpp -lgmp -lmpfr

.PHONY: setup-dependencies test
test: test_main
	./test_main

setup-dependencies:
	curl -L -O https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1.tar.xz
	tar xf CGAL-5.6.1.tar.xz
