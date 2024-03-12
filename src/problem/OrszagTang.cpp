#include "../data.hpp"
#include "../driver.hpp"

void InputParameters::SetParameters()
{
    nx1_   = 256;
    x1min_ = 0.0;
    x1max_ = 2.0 * M_PI;

    nx2_   = 256;
    x2min_ = 0.0;
    x2max_ = 2.0 * M_PI;

    nx3_   = 1;
    x3min_ = 0.0;
    x3max_ = 0.0;

    tlim_       = M_PI;
    delta_tout_ = 1.0e-1;

    dir_name_   = "output/OrszagTang3";

    // recon_type_  = ReconstructionType::DnonerCellSheme;
    // integ_order_ = IntegratorOredr::FirstOrder;

    recon_type_  = ReconstructionType::MusclMinmodScheme;
    integ_order_ = IntegratorOredr::SecondOrder;

    is_set_params_ = true;
}


void InitialCondition(Data *pdata)
{
    int is = pdata->is_, ie = pdata->ie_,
        js = pdata->js_, je = pdata->je_,
        ks = pdata->ks_, ke = pdata->ke_;

    array::Double1D x1c = pdata->x1c_;
    array::Double1D x2c = pdata->x2c_;
    array::Double4D q   = pdata->q_;

    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

                q[IRHO][k][j][i] = GAMMA * GAMMA;
                q[IV1][k][j][i]  = -std::sin(x2c[j]);
                q[IV2][k][j][i]  = std::sin(x1c[i]);
                q[IV3][k][j][i]  = 0.0;
                q[IPR][k][j][i]  = GAMMA;
                q[IB1][k][j][i]  = -std::sin(x2c[j]);
                q[IB2][k][j][i]  = std::sin(2.0*x1c[i]);
                q[IB3][k][j][i]  = 0.0;
                q[IPSI][k][j][i] = 0.0;

            }
        }
    }

    return;
}


void BoundaryCondition(Data *pdata, array::Double4D q)
{
    int is = pdata->is_, ie = pdata->ie_,
        js = pdata->js_, je = pdata->je_,
        ks = pdata->ks_, ke = pdata->ke_;
    int ngh = pdata->ngh_;

    // x1-direction
    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = 1; i <= ngh; ++i) {
                    q[n][k][j][i-1]  = q[n][k][j][ie+i-ngh];
                    q[n][k][j][ie+i] = q[n][k][j][is+i-1];
                }
            }
        }
    }

    // x2-direction
    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int i = is; i <= ie; ++i) {
                for (int j = 1; j <= ngh; ++j) {
                    q[n][k][j-1][i]  = q[n][k][je+j-ngh][i];
                    q[n][k][je+j][i] = q[n][k][js+j-1][i];
                }
            }
        }
    }

    return;
}