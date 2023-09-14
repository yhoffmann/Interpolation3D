.PHONY: test

test:
	g++ -o test -Wall -O3 testing/main.cpp src/*.cpp easy-progress-monitor/src/*.cpp -fopenmp -lgsl