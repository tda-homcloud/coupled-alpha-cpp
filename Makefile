HEADERS := $(wildcard *.hpp)
FLAGS := -Wall -g -std=c++17 -Wall -Ieigen -ICGAL-5.6.1/include
test_simplex: test_simplex.cpp $(HEADERS)
	g++ $(FLAGS) -o test_simplex test_simplex.cpp -lgmp -lmpfr

.PHONY: run_test_simplex
run_test_simplex: test_simplex
	./test_simplex

test_lifted_delaunay: test_lifted_delaunay.cpp $(HEADERS)
	g++ $(FLAGS) -o test_lifted_delaunay test_lifted_delaunay.cpp -lgmp -lmpfr

.PHONY: run_test_lifted_delaunay
run_test_lifted_delaunay: test_lifted_delaunay
	./test_lifted_delaunay

test_min_circumsphere: test_min_circumsphere.cpp $(HEADERS)
	g++ $(FLAGS) -o test_min_circumsphere test_min_circumsphere.cpp -lgmp -lmpfr

.PHONY: run_test_min_circumsphere
run_test_min_circumsphere: test_min_circumsphere
	./test_min_circumsphere

test_relaxed_filtration_value: test_relaxed_filtration_value.cpp $(HEADERS)
	g++ $(FLAGS) -o test_relaxed_filtration_value test_relaxed_filtration_value.cpp -lgmp -lmpfr

.PHONY: run_test_relaxed_filtration_value
run_test_relaxed_filtration_value: test_relaxed_filtration_value
	./test_relaxed_filtration_value

.PHONY: run_test_all
run_test_all:
	make run_test_simplex
	make run_test_lifted_delaunay
	make run_test_min_circumsphere
	make run_test_relaxed_filtration_value

coupled_alpha: main.cpp $(HEADERS)
	g++ $(FLAGS) -o coupled_alpha main.cpp -lgmp -lmpfr

.PHONY: run_cmp
run_cmp: coupled_alpha
	cp compare_with_original_ca.py cmp_work
	cd cmp_work && python compare_with_original_ca.py 3 20 20
	cd cmp_work && python compare_with_original_ca.py 2 20 20

setup-dependencies:
	curl -L -O https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1.tar.xz
	tar xf CGAL-5.6.1.tar.xz
	git clone https://github.com/yohaireani/cycle-registration-persistent-homology.git
	mkdir -p cmp_work
	cp cycle-registration-persistent-homology/*.py cmp_work

