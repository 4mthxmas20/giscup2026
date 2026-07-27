CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pthread

solver: solver.cpp geometry.hpp grid.hpp visibility.hpp
	$(CXX) $(CXXFLAGS) -o $@ solver.cpp

clean:
	rm -f solver

.PHONY: clean
