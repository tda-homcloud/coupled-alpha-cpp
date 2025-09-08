test_main: test_main.cpp coupled_alpha.hpp
	g++ -std=c++17 -Wall -Ieigen -ICGAL-5.6.1/include -o test_main test_main.cpp -lgmp -lmpfr

.PHONY: setup-dependencies
setup-dependencies:
	curl -L -O https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1.tar.xz
	tar xf CGAL-5.6.1.tar.xz
