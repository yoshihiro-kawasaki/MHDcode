#include "data.hpp"
#include "driver.hpp"
#include "timer.hpp"

#include <iostream>

int main(const int argc, const char* argv[])
{
    Timer timer;
    timer.Start();

    InputParameters input;
    input.SetParameters();
    if (! input.is_set_params_) {
        std::cout << "Not set GridInfo" << std::endl;
        return 0;
    }

    Data data(input);

    Driver driver(&data);
    driver.EnrollInitialConditionFunction(InitialCondition);
    driver.EnrollBoundaryConditionFunction(BoundaryCondition);
    driver.SetInitialCondition();

    data.Output(0);

    return 0;

    // main loop 
    long int count = 0;
    double dt, tout;
    
    while (data.time_.t_ < data.time_.tlim_)
    {
        count++;

        // determine dt
        dt = driver.CalculateTimeStep();
        if (data.time_.t_ + dt > data.time_.tlim_) dt = data.time_.tlim_ - data.time_.t_;
        if (data.time_.t_ + dt > tout) dt = tout - data.time_.t_;

        // update
        driver.Integration(dt);
        data.time_.t_ += dt;

        // ouput
        if (data.time_.t_ >= tout || data.time_.t_ == data.time_.tlim_) {
            data.Output(count);
            tout += data.time_.delta_tout_;
        }

    }
    
    timer.End();
    // timer.ShowTimeSeconds();
    timer.ShowTimeMilliSeconds();
    return 0;
}