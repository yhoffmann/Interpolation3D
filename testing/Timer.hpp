#pragma omp


#include <string>
#include <chrono>
#include <iostream>


class Timer
{
    std::string m_scopename;
    std::chrono::time_point<std::chrono::high_resolution_clock,std::chrono::nanoseconds> start;
    std::chrono::time_point<std::chrono::high_resolution_clock,std::chrono::nanoseconds> end;
    
public:

    Timer()
        : start(std::chrono::high_resolution_clock::now())
    {
    }

    Timer(std::string scopename)
        : m_scopename(scopename), start(std::chrono::high_resolution_clock::now())
    {
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
        std::cout << m_scopename << " took " << runtime.count() << " to run\n";
    }
};