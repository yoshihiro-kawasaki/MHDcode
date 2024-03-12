#ifndef DRIVER_HPP_
#define DRIVER_HPP_

// #include "../array.hpp"
// #include "../defs.hpp"
// #include "../data.hpp"

// #include "../eos/eos.hpp"
// #include "../reconstruction/reconstruction.hpp"

#include "array.hpp"
#include "defs.hpp"
#include "data.hpp"

#include "eos/eos.hpp"
#include "reconstruction/reconstruction.hpp"

using InitialConditionFunction  = void(*)(Data *pdata);
using BoundaryConditionFunction = void(*)(Data *pdata, array::Double4D q);

void InitialCondition(Data *pdata);
void BoundaryCondition(Data *pdata, array::Double4D q);

class Driver
{
    friend class Data;
public:
    Driver(Data *pdata, const InputParameters &input);
    ~Driver();

    void EnrollInitialConditionFunction(InitialConditionFunction ic) {
        InitialCondition_ = ic;
    };

    void EnrollBoundaryConditionFunction(BoundaryConditionFunction bc) {
        BoundaryCondition_ = bc;
    };

    void SetInitialCondition() {
        InitialCondition_(pdata_);
        peos_->ConvertPrimToCons(pdata_->q_, pdata_->u_);
    }

    double CalculateTimeStep() {
        return peos_->CalculateTimeStep();
    };

    void Integration(const double dt) {
        (this->*Integration_)(dt);
    }


private:

    void CalculateFluxes(const array::Double4D q);

    void Integration1dProblemFirstOredr(const double dt);
    void Integration1dProblemSecondOredr(const double dt);

    void Integration2dProblemFirstOredr(const double dt);
    void Integration2dProblemSecondOredr(const double dt);

    Data            *pdata_;
    EquationOfState *peos_;
    Reconstruction  *precon_;

    InitialConditionFunction  InitialCondition_;
    BoundaryConditionFunction BoundaryCondition_;

    using IntegrationFunction = void(Driver::*)(const double dt);
    IntegrationFunction Integration_;

    array::Double4D wu_;  // work array of conservative variavles
    array::Double4D wq_;  // work array of primitive variavles

    array::Double4D flx1_; // numerical flux
    array::Double4D flx2_; // numerical flux
    array::Double4D flx3_; // numerical flux

    array::Double2D ql_;
    array::Double2D qlb_;
    array::Double2D qr_;

    double ch_;
};

#endif /* DRIVER_HPP_ */