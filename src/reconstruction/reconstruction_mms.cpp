#include "reconstruction.hpp"

ReconstructionMusclMinmodScheme::ReconstructionMusclMinmodScheme(Data *pdata)
    : Reconstruction(pdata)
{

}


ReconstructionMusclMinmodScheme::~ReconstructionMusclMinmodScheme()
{

}


void ReconstructionMusclMinmodScheme::CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    double dql, dq;
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            dql = q[n][k][j][i]   - q[n][k][j][i-1];
            dq  = q[n][k][j][i+1] - q[n][k][j][i];
            ql[n][i+1] = q[n][k][j][i] + 0.5 * Minmod(dq, dql);
            qr[n][i]   = q[n][k][j][i] - 0.5 * Minmod(dq, dql); 
        }
    }
    return;
}


void ReconstructionMusclMinmodScheme::CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    double dql, dq;
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            dql = q[n][k][j][i]   - q[n][k][j-1][i];
            dq  = q[n][k][j+1][i] - q[n][k][j][i];
            ql[n][i] = q[n][k][j][i] + 0.5 * Minmod(dq, dql);
            qr[n][i] = q[n][k][j][i] - 0.5 * Minmod(dq, dql); 
        }
    }
    return;
}


void ReconstructionMusclMinmodScheme::CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    double dql, dq;
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            dql = q[n][k][j][i]   - q[n][k][j][i-1];
            dq  = q[n][k][j][i+1] - q[n][k][j][i];
            ql[n][i+1] = q[n][k][j][i] + 0.5 * Minmod(dq, dql);
            qr[n][i]   = q[n][k][j][i] - 0.5 * Minmod(dq, dql); 
        }
    }
    return;
}