#include "driver.hpp"

Driver::Driver(Data *pdata)
    : pdata_(pdata), InitialCondition_(nullptr), BoundaryCondition_(nullptr)
{
    int nx1totc = pdata_->nx1totc_, nx1totb = pdata_->nx1totb_,
        nx2totc = pdata_->nx2totc_, nx2totb = pdata_->nx2totb_,
        nx3totc = pdata_->nx3totc_, nx3totb = pdata_->nx3totb_;

    wq_   = array::Allocate4dArray<double>(NVAR, nx3totc, nx2totc, nx1totc);
    wu_   = array::Allocate4dArray<double>(NVAR, nx3totc, nx2totc, nx1totc);
    flx1_ = array::Allocate4dArray<double>(NVAR, nx3totb, nx2totb, nx1totb);
    flx2_ = array::Allocate4dArray<double>(NVAR, nx3totb, nx2totb, nx1totb);
    flx3_ = array::Allocate4dArray<double>(NVAR, nx3totb, nx2totb, nx1totb);
    ql_   = array::Allocate2dArray<double>(NVAR, nx1totb);
    qlb_  = array::Allocate2dArray<double>(NVAR, nx1totb);
    qr_   = array::Allocate2dArray<double>(NVAR, nx1totb);

    // precon_ = new ReconstructionDonorCellScheme(pdata);
    precon_ = new ReconstructionMusclMinmodScheme(pdata_);
    peos_   = new EquationOfState(pdata_);
    Integration_ = &Driver::IntegrationFirstOredr1DProblem;
}


Driver::~Driver()
{
    array::Delete4dArray<double>(wu_);
    array::Delete4dArray<double>(wq_);
    array::Delete4dArray<double>(flx1_);
    array::Delete4dArray<double>(flx2_);
    array::Delete4dArray<double>(flx3_);
    array::Delete2dArray<double>(ql_);
    array::Delete2dArray<double>(qlb_);
    array::Delete2dArray<double>(qr_);

    // delete[] precon_;
    // delete[] peos_;
}


void Driver::CalculateFluxes(const array::Double4D q)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;
    int il, iu, jl, ju, kl, ku; // l:lower , u:upper
    int idir;

     double ch;

    // i-direction
    idir = 1;
    jl   = js;
    ju   = je;
    kl   = ks;
    ku   = ke;
    if (pdata_->nx2_ > 1 && pdata_->nx3_ == 1) {
        jl = js-1;
        ju = je+1;
    } else if (pdata_->nx2_ > 1 && pdata_->nx3_ > 1) {
        jl = js-1;
        ju = je+1;
        kl = ks-1;
        ku = ke+1;
    }

    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            // reconstruct L/R state
            precon_->CalculateReconstructionX1(k, j, is-1, ie+1, q, ql_, qr_);
            // calculate fluxes
            peos_->RiemannSolver(k, j, is, ie+1, idir, ql_, qr_, flx1_, ch);
        }
    }

    // j-direction
    if (pdata_->nx2_ > 1) {

        idir = 2;
        il = is-1;
        iu = ie+1;
        kl = ks;
        ku = ke;

        for (int k = kl; k <= ku; ++k) {

            precon_->CalculateReconstructionX2(k, js-1, il, iu, q, ql_, qr_);

            for (int j = js; j <= je+1; ++j) {

                precon_->CalculateReconstructionX2(k, j, il, iu, q, qlb_, qr_);
                peos_->RiemannSolver(k, j, il, iu, idir, ql_, qr_, flx2_, ch);

                for (int n = 0; n < NVAR; ++n) {
                    for (int i = il; i <= iu; ++i) {
                        ql_[n][i] = qlb_[n][i];
                    }
                }

            }
        }

    }

    // // k-direction
    if (pdata_->nx3_ > 1) {

    }
}


void Driver::IntegrationFirstOredr1DProblem(const double dt)
{
    return;
}



void Driver::IntegrationSecondOredr1DProblem(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, inv_dx1, dudt;

    inv_dx1 = dt / pdata_->dx1_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] - dt * inv_dx1 * (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1;
                    pdata_->u_[n][k][j][i] = 0.5 * pdata_->u_[n][k][j][i] + 0.5 * (wu_[n][k][j][i] + dudt*dt);
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}