#ifndef RECONSTRUCTION_HPP_
#define RECONSTRUCTION_HPP_

#include "../array.hpp"
#include "../defs.hpp"
#include "../data.hpp"

#include <iostream>
#include <cmath>

// class Data;


// #######################################################################
// #######################################################################
// #######################################################################

class Reconstruction
{
public:
    Reconstruction(Data *pdata) : pdata_(pdata) { };
    ~Reconstruction() { };

    virtual void CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr) = 0;
    virtual void CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr) = 0;
    virtual void CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr) = 0;

private:
    
    Data *pdata_;

};

// #######################################################################
// #######################################################################
// #######################################################################


class ReconstructionDonorCellScheme
    : public Reconstruction
{
public:
    ReconstructionDonorCellScheme(Data *pdata);
    ~ReconstructionDonorCellScheme();

    virtual void CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
    virtual void CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
    virtual void CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
private:
};


// #######################################################################
// #######################################################################
// #######################################################################

class ReconstructionMusclMinmodScheme
    : public Reconstruction
{
public:
    ReconstructionMusclMinmodScheme(Data *pdata);
    ~ReconstructionMusclMinmodScheme();

    virtual void CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
    virtual void CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
    virtual void CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr);
private:
    inline double Minmod(double a, double b) {
        return 0.5 * (SIGN(a) + SIGN(b)) * std::min(std::abs(a), std::abs(b));
    }
};


// #######################################################################
// #######################################################################
// #######################################################################

#endif /* RECONSTRUCTION_HPP_ */