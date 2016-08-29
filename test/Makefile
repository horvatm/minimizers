CXX=g++
FC=gfortran

CXXFLAGS=-I../src -Wall -Wextra -O3 -std=c++11

all: test_diff test_minfit test_bfgs_diff test_powell test_praxis test_amoeba

test_diff:
	$(CXX) test_diff.cpp -o test_diff $(CXXFLAGS) -fopenmp

test_minfit:
	$(FC) -c minfit.f90  
	$(CXX) test_minfit.cpp -o test_minfit $(CXXFLAGS) minfit.o -lgfortran 

test_bfgs_diff:
	$(CXX) test_bfgs_diff.cpp -o test_bfgs_diff $(CXXFLAGS) -fopenmp

test_powell:
	$(CXX) test_powell.cpp -o test_powell $(CXXFLAGS)

test_praxis:
	$(CXX) test_praxis.cpp -o test_praxis $(CXXFLAGS)

test_amoeba:
	$(CXX) test_amoeba.cpp -o test_amoeba $(CXXFLAGS)

clean:
	rm -f minfit.o test_diff test_minfit test_bfgs_diff test_powell test_praxis test_amoeba
