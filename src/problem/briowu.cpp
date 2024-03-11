#include "../data.hpp"
#include "../driver.hpp"

void InputParameters::SetParameters()
{
    nx1_   = 64;
    x1min_ = -0.5;
    x1max_ = 0.5;

    nx2_   = 1;
    x2min_ = 0.0;
    x2max_ = 0.0;

    nx3_   = 1;
    x3min_ = 0.0;
    x3max_ = 0.0;

    ngh_   = 2;

    tlim_ = 0.1;

    dir_name_   = "output/briowu";

    is_set_params_ = true;
}


void InitialCondition(Data *pdata)
{
    int is = pdata->is_, ie = pdata->ie_,
        js = pdata->js_, je = pdata->je_,
        ks = pdata->ks_, ke = pdata->ke_;

    array::Double1D x1c = pdata->x1c_;
    array::Double4D q   = pdata->q_;

    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

                if (x1c[i] < 0.0) {
                    q[IRHO][k][j][i] = 1.0;
                    q[IV1][k][j][i]  = 0.0;
                    q[IV2][k][j][i]  = 0.0;
                    q[IV3][k][j][i]  = 0.0;
                    q[IPR][k][j][i]  = 1.0;
                    q[IB1][k][j][i]  = 0.75;
                    q[IB2][k][j][i]  = 1.0;
                    q[IB3][k][j][i]  = 0.0;
                    q[IPSI][k][j][i] = 0.0;
                } else {
                    q[IRHO][k][j][i] = 0.125;
                    q[IV1][k][j][i]  = 0.0;
                    q[IV2][k][j][i]  = 0.0;
                    q[IV3][k][j][i]  = 0.0;
                    q[IPR][k][j][i]  = 0.1;
                    q[IB1][k][j][i]  = 0.75;
                    q[IB2][k][j][i]  = -1.0;
                    q[IB3][k][j][i]  = 0.0;
                    q[IPSI][k][j][i] = 0.0;
                }

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

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = 1; i <= ngh; ++i) {
                    q[n][k][j][is-i] = q[n][k][j][is+i-1];
                }
            }
        }
    }

    for (int n = 0; n < NVAR; ++n) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = 1; i <= ngh; ++i) {
                    q[n][k][j][ie+i] = q[n][k][j][ie-i+1];
                }
            }
        }
    }

    return;
}
