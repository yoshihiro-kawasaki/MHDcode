#include "reconstruction.hpp"

ReconstructionDonorCellScheme::ReconstructionDonorCellScheme(Data *pdata)
    : Reconstruction(pdata)
{

}


ReconstructionDonorCellScheme::~ReconstructionDonorCellScheme()
{

}


void ReconstructionDonorCellScheme::CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            ql[n][i+1] = q[n][k][j][i];
            qr[n][i]   = q[n][k][j][i];
        }
    }
    return;
}


void ReconstructionDonorCellScheme::CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            ql[n][i] = q[n][k][j][i];
            qr[n][i] = q[n][k][j][i];
        }
    }
    return;
}


void ReconstructionDonorCellScheme::CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            ql[n][i] = q[n][k][j][i];
            qr[n][i] = q[n][k][j][i];
        }
    }
    return;
}