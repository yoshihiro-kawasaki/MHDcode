#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <chrono>
#include <iostream>

class Timer
{
public:
    Timer() {

    }

    void Start() {
        tstart_ = std::chrono::system_clock::now();
    }

    void End() {
        tend_ = std::chrono::system_clock::now();
    }

    void ShowTimeMilliSeconds() {
        const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tend_ - tstart_).count());
        std::cout << std::scientific;
        std::cout << "time = " << time << " mil sec" << std::endl;
    }

    void ShowTimeSeconds() {
        const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(tend_ - tstart_).count());
        std::cout << std::scientific;
        std::cout << "time = " << time << " sec" << std::endl;
    }

    void ShowTimeMinutes() {
        const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::minutes>(tend_ - tstart_).count());
        std::cout << std::scientific;
        std::cout << "time = " << time << " min" << std::endl;
    }


    void ShowTimeHours() {
        const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::hours>(tend_ - tstart_).count());
        std::cout << std::scientific;
        std::cout << "time = " << time << " hours" << std::endl;
    }


private:
    
    std::chrono::system_clock::time_point tstart_, tend_;
};

#endif /* TIMER_HPP_ */