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
    ReconstructionType recon_type_ = ReconstructionType::MusclMinmodScheme;
    IntegratorOredr integ_order_ = IntegratorOredr::SecondOrder;
    void SetParameters();
};


class Data
{
public:
    Data(const InputParameters &input);
    ~Data();

    array::Double4D q_;   // primitive    variables (rho, v, p, B)
    array::Double4D u_;   // conservative variables (rho, rhov, e, B)

    array::Double4D flx1_; // numerical flux
    array::Double4D flx2_; // numerical flux
    array::Double4D flx3_; // numerical flux

    array::Double1D x1c_; // cell central
    array::Double1D x2c_; // cell central
    array::Double1D x3c_; // cell central
    array::Double1D x1b_; // cell boundary
    array::Double1D x2b_; // cell boundary
    array::Double1D x3b_; // cell boundary

    int nx1_;
    int nx2_;
    int nx3_;
    int ngh_;
    int nx1totc_; 
    int nx2totc_; 
    int nx3totc_; 
    int nx1totb_; 
    int nx2totb_; 
    int nx3totb_;
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