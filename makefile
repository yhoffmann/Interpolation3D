.PHONY: test

test:
	g++ -Wall -o test testing/main.cpp src/*.cpp -fopenmp