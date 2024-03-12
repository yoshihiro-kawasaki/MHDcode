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

    Driver driver(&data, input);
    driver.EnrollInitialConditionFunction(InitialCondition);
    driver.EnrollBoundaryConditionFunction(BoundaryCondition);
    driver.SetInitialCondition();

    // main loop 
    long int count = 0;
    long int output_count = 0;
    double dt, tout;
    tout = data.time_.delta_tout_;
    data.Output(output_count);
    
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
        // if (data.time_.t_ == data.time_.tlim_) {
            std::cout << "output : " << output_count << std::endl; 
            data.Output(output_count);
            tout += data.time_.delta_tout_;
            output_count++;
        }

    }
    
    timer.End();
    // timer.ShowTimeSeconds();
    timer.ShowTimeMilliSeconds();
    return 0;
}