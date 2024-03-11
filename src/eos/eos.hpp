#ifndef EOS_HPP_
#define EOS_HPP_

#include "../array.hpp"
#include "../defs.hpp"
#include "../data.hpp"

#include <iostream>
#include <cmath>

class EquationOfState
{
public:
    EquationOfState(Data *pdata);
    ~EquationOfState();

    virtual void ConvertConsToPrim(const array::Double4D u, array::Double4D q);
    virtual void ConvertPrimToCons(const array::Double4D q, array::Double4D u);
    virtual double CalculateTimeStep();
    virtual void RiemannSolver(const int k, const int j, const int is, const int ie, const int idir, 
        array::Double2D ql, array::Double2D qr, array::Double4D flx, double ch);

    // void ConvertConsToPrim(const array::Double4D u, array::Double4D q);
    // void ConvertPrimToCons(const array::Double4D q, array::Double4D u);
    // double CalculateTimeStep();
    // void RiemannSolver(const int k, const int j, const int is, const int ie, const int idir, 
    //     array::Double2D ql, array::Double2D qr, array::Double4D flx, double ch);

private:
    Data *pdata_;
};

#endif /* EOS_HPP_ */