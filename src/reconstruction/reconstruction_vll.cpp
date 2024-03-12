#include "reconstruction.hpp"

ReconstructionMusclVanLeerScheme::ReconstructionMusclVanLeerScheme(Data *pdata)
    : Reconstruction(pdata)
{
    int nx1totv = pdata_->nx1totv_, nx1totf = pdata_->nx1totf_,
        nx2totv = pdata_->nx2totv_, nx2totf = pdata_->nx2totf_,
        nx3totv = pdata_->nx3totv_, nx3totf = pdata_->nx3totf_;

    dql_  = array::Allocate2dArray<double>(NVAR, nx1totf);
    dqr_  = array::Allocate2dArray<double>(NVAR, nx1totf);
    dqm_  = array::Allocate2dArray<double>(NVAR, nx1totf);
    qc_   = array::Allocate2dArray<double>(NVAR, nx1totf);
    dx1p_ = array::Allocate1dArray<double>(nx1totv);
    dx1m_ = array::Allocate1dArray<double>(nx1totv);

    int is = pdata_->is_, ie = pdata_->ie_;
    for (int i = is; i <= ie; ++i) {
        dx1p_[i] = (pdata_->x1f_[i+1] - pdata_->x1v_[i]) / pdata_->dx1f_[i];
        dx1m_[i] = (pdata_->x1v_[i]   - pdata_->x1f_[i]) / pdata_->dx1f_[i];
    }

}


ReconstructionMusclVanLeerScheme::~ReconstructionMusclVanLeerScheme()
{
    array::Delete2dArray<double>(dql_);
    array::Delete2dArray<double>(dqr_);
    array::Delete2dArray<double>(dqm_);
    array::Delete2dArray<double>(qc_);
    array::Delete1dArray<double>(dx1p_);
    array::Delete1dArray<double>(dx1m_);
}


void ReconstructionMusclVanLeerScheme::CalculateReconstructionX1(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            dql_[n][i] = q[n][k][j][i]   - q[n][k][j][i-1];
            dqr_[n][i] = q[n][k][j][i+1] - q[n][k][j][i];
            qc_[n][i]  = q[n][k][j][i];
        }
    }

    double dq2;
    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            dq2 = dql_[n][i] * dqr_[n][i];
            dqm_[n][i] = (dq2 > 0.0 ? 2.0 * dq2 / (dql_[n][i] + dqr_[n][i]) : 0.0);
        }
    }

    for (int n = 0; n < NVAR; ++n) {
        for (int i = il; i <= iu; ++i) {
            ql[n][i+1] = qc_[n][i] + dx1p_[i] * dqm_[n][i];
            qr[n][i]   = qc_[n][i] - dx1m_[i] * dqm_[n][i];
        }
    }
}


void ReconstructionMusclVanLeerScheme::CalculateReconstructionX2(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{

}


void ReconstructionMusclVanLeerScheme::CalculateReconstructionX3(const int k, const int j, const int il, const int iu,
        const array::Double4D q, array::Double2D ql, array::Double2D qr)
{

}