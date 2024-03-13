#include "driver.hpp"

Driver::Driver(Data *pdata, const InputParameters &input)
    : pdata_(pdata), precon_{nullptr, nullptr}, peos_(nullptr), InitialCondition_(nullptr), BoundaryCondition_(nullptr), Integration_(nullptr)
{
    int nx1totv = pdata_->nx1totv_, nx1totf = pdata_->nx1totf_,
        nx2totv = pdata_->nx2totv_, nx2totf = pdata_->nx2totf_,
        nx3totv = pdata_->nx3totv_, nx3totf = pdata_->nx3totf_;

    wq_   = array::Allocate4dArray<double>(NVAR, nx3totv, nx2totv, nx1totv);
    wu_   = array::Allocate4dArray<double>(NVAR, nx3totv, nx2totv, nx1totv);
    flx1_ = array::Allocate4dArray<double>(NVAR, nx3totf, nx2totf, nx1totf);
    flx2_ = array::Allocate4dArray<double>(NVAR, nx3totf, nx2totf, nx1totf);
    flx3_ = array::Allocate4dArray<double>(NVAR, nx3totf, nx2totf, nx1totf);
    ql_   = array::Allocate2dArray<double>(NVAR, nx1totf);
    qlb_  = array::Allocate2dArray<double>(NVAR, nx1totf);
    qr_   = array::Allocate2dArray<double>(NVAR, nx1totf);

    if (input.recon_type_ == ReconstructionType::DnonerCellSheme) {
        precon_[0] = new ReconstructionDonorCellScheme(pdata);
    } else if (input.recon_type_ == ReconstructionType::MinmodScheme) {
        precon_[0] = new ReconstructionMusclMinmodScheme(pdata_);
    } else if (input.recon_type_ == ReconstructionType::VanLeerScheme) {
        precon_[0] = new ReconstructionMusclVanLeerScheme(pdata_);
    }

    peos_   = new EquationOfState(pdata_);

    if (pdata_->xdim_ == DimensionsOfProblem::One) {
        if (input.integ_type_ == IntegratorType::ForwardEuler) {
            Integration_ = &Driver::Integration1dProblemForwardEuler;
        } else if (input.integ_type_ == IntegratorType::SSPRK22) {
            Integration_ = &Driver::Integration1dProblemSSPRK22;
        } else if (input.integ_type_ == IntegratorType::PCM) {
            Integration_ = &Driver::Integration1dProblemPCM;
            precon_[1] = new ReconstructionDonorCellScheme(pdata_);
        }
    } else if (pdata_->xdim_ == DimensionsOfProblem::Two) {
        if (input.integ_type_ == IntegratorType::ForwardEuler) {
            Integration_ = &Driver::Integration2dProblemForwardEuler;
        } else if (input.integ_type_ == IntegratorType::SSPRK22) {
            Integration_ = &Driver::Integration2dProblemSSPRK22;
        } else if (input.integ_type_ == IntegratorType::PCM) {
            // Integration_ = 
        }
    } else if (pdata_->xdim_ == DimensionsOfProblem::Three) {

    }

    // check
    if (precon_ == nullptr) {
        std::cout << "precon == nullptr" << std::endl;
        std::exit(1);
    }

    if (peos_ == nullptr) {
        std::cout << "peos == nullptr" << std::endl;
        std::exit(1);
    }

    if (Integration_ == nullptr) {
        std::cout << "Integration == nullptr" << std::endl;
        std::exit(1);
    }
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
}



void Driver::CalculateFluxes(const array::Double4D q, Reconstruction *precon, EquationOfState *peos)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;
    int il, iu, jl, ju, kl, ku; // l:lower , u:upper
    int idir;
    double ch = ch_;

    // i-direction
    idir = CoordinateDirection::X1dir;
    jl   = js;
    ju   = je;
    kl   = ks;
    ku   = ke;
    if (pdata_->nx2_ > 1 && pdata_->nx3_ == 1) { 
        // 2D
        jl = js-1;
        ju = je+1;
    } else if (pdata_->nx2_ > 1 && pdata_->nx3_ > 1) {
        // 3D
        jl = js-1;
        ju = je+1;
        kl = ks-1;
        ku = ke+1;
    }

    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            // reconstruct L/R state
            precon->CalculateReconstructionX1(k, j, is-1, ie+1, q, ql_, qr_);
            // calculate fluxes
            peos->RiemannSolver(k, j, is, ie+1, idir, ql_, qr_, flx1_, ch);
        }
    }

    // j-direction
    if (pdata_->xdim_ == DimensionsOfProblem::Two) {

        idir = CoordinateDirection::X2dir;
        il = is-1;
        iu = ie+1;
        kl = ks;
        ku = ke;
        if (pdata_->nx3_ > 1) {
            // 3D
            kl = ks - 1;
            ku = kl - 1;
        }

        for (int k = kl; k <= ku; ++k) {

            precon->CalculateReconstructionX2(k, js-1, il, iu, q, ql_, qr_);

            for (int j = js; j <= je+1; ++j) {

                precon->CalculateReconstructionX2(k, j, il, iu, q, qlb_, qr_);
                peos->RiemannSolver(k, j, il, iu, idir, ql_, qr_, flx2_, ch);

                for (int n = 0; n < NVAR; ++n) {
                    for (int i = il; i <= iu; ++i) {
                        ql_[n][i] = qlb_[n][i];
                    }
                }

            }
        }

    }

    // k-direction
    if (pdata_->xdim_ == DimensionsOfProblem::Three) {

        idir = CoordinateDirection::X3dir;
        il = is - 1;
        iu = ie + 1;
        jl = js - 1;
        ju = je + 1;

        for (int j = jl; j <= ju; ++j) {

            precon->CalculateReconstructionX3(ks-1, j, il, iu, q, ql_, qr_);

            for (int k = ks; k <= ke+1; ++k) {

                precon->CalculateReconstructionX3(k, j, il, iu, q, qlb_, qr_);
                peos->RiemannSolver(k, j, il, iu, idir, ql_, qr_, flx3_, ch);

                for (int n = 0; n < NVAR; ++n) {
                    for (int i = il; i <= iu; ++i) {
                        ql_[n][i] = qlb_[n][i];
                    }
                }

            }
        }
    }

    return;
}


void Driver::Integration1dProblemForwardEuler(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, inv_dx1, dudt;
    ch_ = CFL * pdata_->dx1_ / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i] - dt * inv_dx1 * (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]);
                    if (n == IPSI) pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);
    
    return;
}



void Driver::Integration1dProblemSSPRK22(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, inv_dx1, dudt;
    ch_ = CFL * pdata_->dx1_ / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

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
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1;
                    pdata_->u_[n][k][j][i] = 0.5 * (pdata_->u_[n][k][j][i] + wu_[n][k][j][i] + dudt*dt);
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}



void Driver::Integration1dProblemPCM(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, inv_dx1, dudt;
    ch_ = CFL * pdata_->dx1_ / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[1], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] - 0.5 * dt * inv_dx1 * (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1;
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i] + dudt*dt;
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}


void Driver::Integration2dProblemForwardEuler(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx2, inv_dx1, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                            - dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                                  + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2);
                    if (n == IPSI) pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);
    
    return;
}



void Driver::Integration2dProblemSSPRK22(const double dt)
{
int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx2, inv_dx1, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                    - dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                          + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                           - (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2;
                    pdata_->u_[n][k][j][i] = 0.5 * pdata_->u_[n][k][j][i] + 0.5 * (wu_[n][k][j][i] + dudt*dt);
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}


void Driver::Integration2dProblemPCM(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx2, inv_dx1, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[1], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                    - 0.5 * dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                                + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                           - (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2;
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i]  + dudt*dt;
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}


void Driver::Integration3dProblemForwardEuler(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx1, inv_dx2, inv_dx3, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    dxmin = std::min(pdata_->dx3_, dxmin);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;
    inv_dx3 = 1.0 / pdata_->dx3_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                            - dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                                  + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2
                                                  + (flx3_[n][k+1][j][i] - flx3_[n][k][j][i]) * inv_dx3);
                    if (n == IPSI) pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);
    
    return;
}



void Driver::Integration3dProblemSSPRK22(const double dt)
{
int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx1, inv_dx2, inv_dx3, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    dxmin = std::min(pdata_->dx3_, dxmin);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;
    inv_dx3 = 1.0 / pdata_->dx3_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                    - dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                          + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2
                                          + (flx3_[n][k+1][j][i] - flx3_[n][k][j][i]) * inv_dx3);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                           - (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2
                           - (flx3_[n][k+1][j][i] - flx3_[n][k][j][i]) * inv_dx3;
                    pdata_->u_[n][k][j][i] = 0.5 * pdata_->u_[n][k][j][i] + 0.5 * (wu_[n][k][j][i] + dudt*dt);
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}


void Driver::Integration3dProblemPCM(const double dt)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double cd, dxmin, inv_dx1, inv_dx2, inv_dx3, dudt;

    dxmin = std::min(pdata_->dx1_, pdata_->dx2_);
    dxmin = std::min(pdata_->dx3_, dxmin);
    ch_ = CFL * dxmin / dt;
    cd  = std::exp(-dt*ch_/CR);
    inv_dx1 = 1.0 / pdata_->dx1_;
    inv_dx2 = 1.0 / pdata_->dx2_;
    inv_dx3 = 1.0 / pdata_->dx3_;

    BoundaryCondition_(pdata_, pdata_->q_);
    CalculateFluxes(pdata_->q_, precon_[1], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    wu_[n][k][j][i] = pdata_->u_[n][k][j][i] 
                                    - 0.5 * dt * ((flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                                                + (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2
                                                + (flx3_[n][k+1][j][i] - flx3_[n][k][j][i]) * inv_dx3);
                    if (n == IPSI) wu_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(wu_, wq_);
    BoundaryCondition_(pdata_, wq_);
    CalculateFluxes(wq_, precon_[0], peos_);

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dudt = - (flx1_[n][k][j][i+1] - flx1_[n][k][j][i]) * inv_dx1
                           - (flx2_[n][k][j+1][i] - flx2_[n][k][j][i]) * inv_dx2
                           - (flx3_[n][k+1][j][i] - flx3_[n][k][j][i]) * inv_dx3;
                    pdata_->u_[n][k][j][i] = pdata_->u_[n][k][j][i]  + dudt*dt;
                    if (n == IPSI)  pdata_->u_[n][k][j][i] *= cd;
                }
            }
        }
    }

    peos_->ConvertConsToPrim(pdata_->u_, pdata_->q_);

    return;
}


