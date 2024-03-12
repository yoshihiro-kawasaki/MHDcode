#ifndef DATA_HPP_
#define DATA_HPP_

#include "array.hpp"
#include "defs.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>


struct InputParameters
{
    double x1min_;
    double x2min_;
    double x3min_;
    double x1max_;
    double x2max_;
    double x3max_;
    int    nx1_;
    int    nx2_;
    int    nx3_;
    double tlim_;
    double delta_tout_;
    std::string dir_name_;
    bool   is_set_params_ = false;
    ReconstructionType recon_type_ = ReconstructionType::MinmodScheme;
    IntegratorType integ_type_ = IntegratorType::SSPRK22;
    void SetParameters();
};


class Data
{
public:
    Data(const InputParameters &input);
    ~Data();

    array::Double4D q_;   // primitive    variables (rho, v, p, B)
    array::Double4D u_;   // conservative variables (rho, rhov, e, B)

    array::Double1D x1v_;  // cell volume averaged position
    array::Double1D x2v_;  // cell volume averaged position
    array::Double1D x3v_;  // cell volume averaged position
    array::Double1D dx1v_; // cell volume averaged spacing
    array::Double1D dx2v_; // cell volume averaged spacing
    array::Double1D dx3v_; // cell volume averaged spacing
    array::Double1D x1f_;  // cell face position
    array::Double1D x2f_;  // cell face position
    array::Double1D x3f_;  // cell face position
    array::Double1D dx1f_; // cell face spacing
    array::Double1D dx2f_; // cell face sapcing
    array::Double1D dx3f_; // cell face spacing

    int nx1_;
    int nx2_;
    int nx3_;
    int ngh_;
    int nx1totv_; 
    int nx2totv_; 
    int nx3totv_; 
    int nx1totf_; 
    int nx2totf_; 
    int nx3totf_;
    double x1min_;
    double x2min_;
    double x3min_;
    double x1max_;
    double x2max_;
    double x3max_;
    double dx1_;
    double dx2_;
    double dx3_;
    DimensionsOfProblem xdim_ = DimensionsOfProblem::One;

    int is_;
    int js_;
    int ks_;
    int ie_;
    int je_;
    int ke_;

    struct Time
    {
        double t_;
        double tlim_;
        double delta_tout_;
    } time_;

    std::string dir_name_;

    void Output(const long int count);
    void OutputGrid();

protected:

private:

    bool is_allocate_array_;
};

#endif /* DATA_HPP_ */