#include "eos.hpp"

EquationOfState::EquationOfState(Data *pdata)
    : pdata_(pdata)
{

}


EquationOfState::~EquationOfState()
{

}


void EquationOfState::ConvertConsToPrim(const array::Double4D u, array::Double4D q)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double gam_m1 = GAMMA - 1.0;
    double inv_rho;
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

                q[IRHO][k][j][i] = u[IRHO][k][j][i];
                inv_rho          = 1.0 / q[IRHO][k][j][i];
                q[IV1][k][j][i]  = u[IM1][k][j][i] * inv_rho;
                q[IV2][k][j][i]  = u[IM2][k][j][i] * inv_rho;
                q[IV3][k][j][i]  = u[IM3][k][j][i] * inv_rho;

                q[IB1][k][j][i]  = u[IB1][k][j][i];
                q[IB2][k][j][i]  = u[IB2][k][j][i];
                q[IB3][k][j][i]  = u[IB3][k][j][i];

                q[IPR][k][j][i]  = gam_m1 * (u[IEN][k][j][i]
                                    - 0.5 * (q[IRHO][k][j][i] * (SQR(q[IV1][k][j][i]) + SQR(q[IV2][k][j][i]) + SQR(q[IV3][k][j][i])) 
                                        + SQR(q[IB1][k][j][i]) + SQR(q[IB2][k][j][i]) + SQR(q[IB3][k][j][i])));

                q[IPSI][k][j][i] = u[IPSI][k][j][i];

            }
        }
    }
}



void EquationOfState::ConvertPrimToCons(const array::Double4D q, array::Double4D u)
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double inv_gam_m1 = GAMMA - 1.0;

    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

                u[IRHO][k][j][i] = q[IRHO][k][j][i];
                u[IM1][k][j][i]  = q[IRHO][k][j][i] * q[IV1][k][j][i];
                u[IM2][k][j][i]  = q[IRHO][k][j][i] * q[IV2][k][j][i];
                u[IM3][k][j][i]  = q[IRHO][k][j][i] * q[IV3][k][j][i];

                u[IEN][k][j][i]  = 0.5 * q[IRHO][k][j][i] * (SQR(q[IV1][k][j][i]) + SQR(q[IV2][k][j][i]) + SQR(q[IV3][k][j][i]))
                                    + 0.5 * (SQR(q[IB1][k][j][i]) + SQR(q[IB2][k][j][i]) + SQR(q[IB3][k][j][i]))
                                    + q[IPR][k][j][i] * inv_gam_m1;

                u[IB1][k][j][i]  = q[IB1][k][j][i];
                u[IB2][k][j][i]  = q[IB2][k][j][i];
                u[IB3][k][j][i]  = q[IB3][k][j][i];

                u[IPSI][k][j][i] = q[IPSI][k][j][i];

            }
        }
    }
}

double EquationOfState::CalculateTimeStep()
{
    int is = pdata_->is_, ie = pdata_->ie_,
        js = pdata_->js_, je = pdata_->je_,
        ks = pdata_->ks_, ke = pdata_->ke_;

    double dt, dx1, dx2, dx3;
    double rho, inv_rho, v1, v2, v3, pr, b1, b2, b3;
    double cssq, vasq, vf1, vf2, vf3;
    array::Double4D q = pdata_->q_;
    dt  = 1e300;
    dx1 = pdata_->dx1_;
    dx2 = pdata_->dx2_;
    dx3 = pdata_->dx3_;

    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {

                rho = q[IRHO][k][j][i];
                v1  = q[IV1][k][j][i];
                v2  = q[IV2][k][j][i];
                v3  = q[IV3][k][j][i];
                pr  = q[IPR][k][j][i];
                b1  = q[IB1][k][j][i];
                b2  = q[IB2][k][j][i];
                b3  = q[IB3][k][j][i];

                inv_rho = 1.0 / rho;
                cssq    = GAMMA * pr * inv_rho;
                vasq    = (b1*b1 + b2*b2 + b3*b3) * inv_rho;

                // calculate fast magnetosinic speed
                vf1  = std::sqrt(0.5 * (cssq + vasq + std::sqrt(SQR(cssq - vasq) + 4.0*cssq*(b2*b2 + b3*b3)*inv_rho)));
                vf2  = std::sqrt(0.5 * (cssq + vasq + std::sqrt(SQR(vasq - cssq) + 4.0*cssq*(b3*b3 + b1*b1)*inv_rho)));
                vf3  = std::sqrt(0.5 * (cssq + vasq + std::sqrt(SQR(vasq - cssq) + 4.0*cssq*(b1*b1 + b2*b2)*inv_rho)));

                dt   = std::min(dt, dx1/(std::abs(v1) + vf1));
                dt   = (pdata_->nx2_ > 1 ? std::min(dt, dx2/(std::abs(v2) + vf2)) : dt);
                dt   = (pdata_->nx3_ > 1 ? std::min(dt, dx3/(std::abs(v3) + vf3)) : dt);
            }
        }
    }

    return CFL*dt;
}

/*
    Hlld RiemannSolver
*/
void EquationOfState::RiemannSolver(const int k, const int j, const int il, const int iu, const int idir, 
        array::Double2D ql, array::Double2D qr, array::Double4D flx, double ch)
{
    int IVN, IVT, IVU, IBN, IBT, IBU;
    double wl[NVAR], wr[NVAR];
    double ul[NVAR], ur[NVAR];
    double ulst[NVAR], urst[NVAR];
    double uldst[NVAR], urdst[NVAR];
    double fl[NVAR], fr[NVAR], fw[NVAR];

    double gam        = GAMMA;
    double inv_gam_m1 = 1.0 / (gam - 1.0);

    if (idir == 1) {
        IVN = IV1;
        IVT = IV2;
        IVU = IV3;
        IBN = IB1;
        IBT = IB2;
        IBU = IB3;
    } else if (idir == 2) {
        IVN = IV2;
        IVT = IV3;
        IVU = IV1;
        IBN = IB2;
        IBT = IB3;
        IBU = IB1;
    } else if (idir == 3) {
        IVN = IV3;
        IVT = IV1;
        IVU = IV2;
        IBN = IB3;
        IBT = IB1;
        IBU = IB2;
    }

    for (int i = il; i <= iu; ++i) {

        // step 1
        wl[IRHO] = ql[IRHO][i];
        wl[IV1]  = ql[IVN][i];
        wl[IV2]  = ql[IVT][i];
        wl[IV3]  = ql[IVU][i];
        wl[IPR]  = ql[IPR][i];
        wl[IB1]  = ql[IBN][i];
        wl[IB2]  = ql[IBT][i];
        wl[IB3]  = ql[IBU][i];
        wl[IPSI] = ql[IPSI][i];

        wr[IRHO] = qr[IRHO][i];
        wr[IV1]  = qr[IVN][i];
        wr[IV2]  = qr[IVT][i];
        wr[IV3]  = qr[IVU][i];
        wr[IPR]  = qr[IPR][i];
        wr[IB1]  = qr[IBN][i];
        wr[IB2]  = qr[IBT][i];
        wr[IB3]  = qr[IBU][i];
        wr[IPSI] = qr[IPSI][i];

        double inv_rhol   = 1.0 / wl[IRHO];
        double inv_rhor   = 1.0 / wr[IRHO];

        // 9-wave (D2002 (42))
        double bxm  = wl[IB1]  + 0.5*(wr[IB1]  - wl[IB1])  - 0.5*(wr[IPSI] - wl[IPSI])/ch;
        double psim = wl[IPSI] + 0.5*(wr[IPSI] - wl[IPSI]) - 0.5*ch*(wr[IB1] - wl[IB1]);

        double bx_sq  = bxm*bxm;
        double abs_bx = std::abs(bxm);

        // magnetic pressure
        double pbl = 0.5 * (bx_sq + wl[IB2]*wl[IB2] + wl[IB3]*wl[IB3]);
        double pbr = 0.5 * (bx_sq + wr[IB2]*wr[IB2] + wr[IB3]*wr[IB3]);
        // total pressure (= thermal pressure + magnetic pressure)
        double ptl = wl[IPR] + pbl;
        double ptr = wr[IPR] + pbr;

        // conserved variables
        ul[IRHO] = wl[IRHO];
        ul[IM1]  = wl[IRHO]*wl[IV1];
        ul[IM2]  = wl[IRHO]*wl[IV2];
        ul[IM3]  = wl[IRHO]*wl[IV3];
        ul[IEN]  = inv_gam_m1*wl[IPR] + pbl + 0.5*wl[IRHO]*(wl[IV1]*wl[IV1] + wl[IV2]*wl[IV2] + wl[IV3]*wl[IV3]);
        ul[IB2]  = wl[IB2];
        ul[IB3]  = wl[IB3];

        ur[IRHO] = wr[IRHO];
        ur[IM1]  = wr[IRHO]*wr[IV1];
        ur[IM2]  = wr[IRHO]*wr[IV2];
        ur[IM3]  = wr[IRHO]*wr[IV3];
        ur[IEN]  = inv_gam_m1*wr[IPR] + pbr + 0.5*wr[IRHO]*(wr[IV1]*wr[IV1] + wr[IV2]*wr[IV2] + wr[IV3]*wr[IV3]);
        ur[IB2]  = wr[IB2];
        ur[IB3]  = wr[IB3];

        // step 2
        // magneto-acoustic wave speed (fast-mode)
        double gamprl = gam * wl[IPR];
        double gamprr = gam * wr[IPR];
        double cfl = std::sqrt(0.5*inv_rhol*( (gamprl + 2.0*pbl) + std::sqrt( (gamprl - 2.0*pbl)*(gamprl - 2.0*pbl) 
            + 4.0*gamprl*( wl[IB2]*wl[IB2] + wl[IB3]*wl[IB3] ) ) ) );
        double cfr = std::sqrt( 0.5*inv_rhor*( (gamprr + 2.0*pbr) + std::sqrt( (gamprr - 2.0*pbr)*(gamprr - 2.0*pbr) 
            + 4.0*gamprr*( wr[IB2]*wr[IB2] + wr[IB3]*wr[IB3] ) ) ) );

        double s0 = std::min( wl[IV1] - cfl, wr[IV1] - cfr );
        double s4 = std::max( wl[IV1] + cfl, wr[IV1] + cfr );

        fl[IRHO] = ul[IM1];
        fl[IM1]  = ul[IM1] * wl[IV1] + ptl - bx_sq;
        fl[IM2]  = ul[IM2] * wl[IV1] - bxm * wl[IB2];
        fl[IM3]  = ul[IM3] * wl[IV1] - bxm * wl[IB3];
        fl[IEN]  = wl[IV1] * (ul[IEN] + ptl - bx_sq) - bxm * (wl[IV2] * ul[IB2] + wl[IV3] * ul[IB3]);
        fl[IB2]  = ul[IB2] * wl[IV1] - bxm * wl[IV2];
        fl[IB3]  = ul[IB3] * wl[IV1] - bxm * wl[IV3];

        fr[IRHO] = ur[IM1];
        fr[IM1]  = ur[IM1] * wr[IV1] + ptr - bx_sq;
        fr[IM2]  = ur[IM2] * wr[IV1] - bxm * wr[IB2];
        fr[IM3]  = ur[IM3] * wr[IV1] - bxm * wr[IB3];
        fr[IEN]  = wr[IV1] * (ur[IEN] + ptr - bx_sq) - bxm * (wr[IV2] * ur[IB2] + wr[IV3] * ur[IB3]);
        fr[IB2]  = ur[IB2] * wr[IV1] - bxm * wr[IV2];
        fr[IB3]  = ur[IB3] * wr[IV1] - bxm * wr[IV3];

        // step 4
        double s0_v1l = s0 - wl[IV1];
        double s4_v1r = s4 - wr[IV1];

        // S_M: eqn (38) of Miyoshi & Kusano
        double s2 = (s4_v1r * ur[IM1] - s0_v1l * ul[IM1] + (ptl - ptr)) / (s4_v1r * ur[IRHO] - s0_v1l * ul[IRHO]);

        double s0_s2 = s0 - s2;
        double s4_s2 = s4 - s2;
        double inv_s0_s2 = 1.0 / s0_s2;
        double inv_s4_s2 = 1.0 / s4_s2;
        // eqn (43) of Miyoshi & Kusano
        ulst[IRHO] = ul[IRHO] * s0_v1l * inv_s0_s2;
        urst[IRHO] = ur[IRHO] * s4_v1r * inv_s4_s2;
        double inv_ulst_d = 1.0 / ulst[IRHO];
        double inv_urst_d = 1.0 / urst[IRHO];
        double sqrt_ulst_d = std::sqrt(ulst[IRHO]);
        double sqrt_urst_d = std::sqrt(urst[IRHO]);

        // eqn (51) of Miyoshi & Kusano
        double s1 = s2 - abs_bx / sqrt_ulst_d;
        double s3 = s2 + abs_bx / sqrt_urst_d;

        // step 5
        // eqn (23) explicitly becomes eq (41) of Miyoshi & Kusano
        double ptstl = ptl + ul[IRHO] * s0_v1l * (s2 - wl[IV1]);
        double ptstr = ptr + ur[IRHO] * s4_v1r * (s2 - wr[IV1]);
        double ptst  = 0.5 * (ptstl + ptstr);

        // ul* - eqn (39) of M&K
        ulst[IM1] = ulst[IRHO] * s2;
        double denoml_st = ul[IRHO] * s0_v1l * s0_s2 - bx_sq;
        if (std::abs(denoml_st) < 1.0e-8*ptst) {
            // Degenerate case
            ulst[IM2] = ulst[IRHO] * wl[IV2];
            ulst[IM3] = ulst[IRHO] * wl[IV3];
            ulst[IB2] = ul[IB2];
            ulst[IB3] = ul[IB3];
        } else {
            denoml_st = 1.0 / denoml_st;

            // eqns (44) and (46) of M&K
            double tmp = bxm*(s0_v1l - s0_s2) * denoml_st;
            ulst[IM2] = ulst[IRHO] * (wl[IV2] - ul[IB2] * tmp);
            ulst[IM3] = ulst[IRHO] * (wl[IV3] - ul[IB3] * tmp);
            // eqns (45) and (47) of M&K
            tmp = (ul[IRHO] * s0_v1l * s0_v1l - bx_sq) * denoml_st;
            ulst[IB2] = ul[IB2] * tmp;
            ulst[IB3] = ul[IB3] * tmp;
        }
        // vlst dot blst
        double vbstl = (ulst[IM1] * bxm + ulst[IM2] * ulst[IB2] + ulst[IM3] * ulst[IB3]) * inv_ulst_d;
        // eqn (48) of M&K
        ulst[IEN] = (s0_v1l * ul[IEN] - ptl * wl[IV1] + ptst*s2 + bxm*(wl[IV1] * bxm + wl[IV2] * ul[IB2]
            + wl[IV3] * ul[IB3] - vbstl)) * inv_s0_s2;

        // ur* - eqn (39) of M&K
        urst[IM1] = urst[IRHO] * s2;
        double denomr_st = ur[IRHO] * s4_v1r * s4_s2 - bx_sq;
        if (std::abs(denomr_st) < 1.0e-8*ptst) {
            // Degenerate case
            urst[IM2] = urst[IRHO] * wr[IV2];
            urst[IM3] = urst[IRHO] * wr[IV3];
            urst[IB2] = ur[IB2];
            urst[IB3] = ur[IB3];
        } else {
            denomr_st = 1.0 / denomr_st;
            // eqns (44) and (46) of M&K
            double tmp = bxm * (s4_v1r - s4_s2) * denomr_st;
            urst[IM2] = urst[IRHO] * (wr[IV2] - ur[IB2] * tmp);
            urst[IM3] = urst[IRHO] * (wr[IV3] - ur[IB3] * tmp);
            // eqns (45) and (47) of M&K
            tmp = (ur[IRHO] * s4_v1r * s4_v1r - bx_sq) * denomr_st;
            urst[IB2] = ur[IB2] * tmp;
            urst[IB3] = ur[IB3] * tmp;
        }
        // vrst dot brst
        double vbstr = (urst[IM1] * bxm + urst[IM2] * urst[IB2] + urst[IM3] * urst[IB3]) * inv_urst_d;
        urst[IEN] = (s4_v1r * ur[IEN] - ptr*wr[IV1] + ptst*s2 + bxm * (wr[IV1] * bxm + wr[IV2] * ur[IB2]
            + wr[IV3] * ur[IB3] - vbstr)) * inv_s4_s2;
            
        // ul** and ur** - if Bx is near zero, same as *-states
        if (0.5*bx_sq < 1.0e-8*ptst) {
            for (int n = 0; n < NVAR; ++n) {
                uldst[n] = ulst[n];
                urdst[n] = urst[n];
            }
        } else {
            double inv_sqrt_sum_ust_d = 1.0 / (sqrt_ulst_d + sqrt_urst_d);
            double sgn_bx = SIGN(bxm);

            uldst[IRHO] = ulst[IRHO];
            urdst[IRHO] = urst[IRHO];

            uldst[IV1] = ulst[IV1];
            urdst[IV1] = urst[IV1];

            // eqn (59) of M&K
            double tmp = inv_sqrt_sum_ust_d * (sqrt_ulst_d * ulst[IM2] * inv_ulst_d + sqrt_urst_d * urst[IM2] * inv_urst_d
                + sgn_bx * (urst[IB2] - ulst[IB2]));
            uldst[IM2] = uldst[IRHO] * tmp;
            urdst[IM2] = urdst[IRHO] * tmp;

            // eqn (60) of M&K
            tmp = inv_sqrt_sum_ust_d * (sqrt_ulst_d * ulst[IM3] * inv_ulst_d + sqrt_urst_d * urst[IM3] * inv_urst_d 
                + sgn_bx * (urst[IB3] - ulst[IB3]));
            uldst[IM3] = uldst[IRHO] * tmp;
            urdst[IM3] = urdst[IRHO] * tmp;

            // eqn (61) of M&K
            tmp = inv_sqrt_sum_ust_d * (sqrt_ulst_d * urst[IB2] + sqrt_urst_d * ulst[IB2] 
                + sgn_bx * sqrt_ulst_d * sqrt_urst_d * (urst[IM2] * inv_urst_d - ulst[IM2] * inv_ulst_d));
            uldst[IB2] = tmp;
            urdst[IB2] = tmp;

            // eqn (62) of M&K
            tmp = inv_sqrt_sum_ust_d * (sqrt_ulst_d * urst[IB3] + sqrt_urst_d * ulst[IB3] 
                + sgn_bx * sqrt_ulst_d * sqrt_urst_d * (urst[IM3] * inv_urst_d - ulst[IM3] * inv_ulst_d));
            uldst[IB3] = tmp;
            urdst[IB3] = tmp;

            // eqn (63) of M&K
            tmp = s2 * bxm + (uldst[IM2] * uldst[IB2] + uldst[IM3] * uldst[IB3]) / uldst[IRHO];
            uldst[IEN] = ulst[IEN] - sqrt_ulst_d * sgn_bx * (vbstl - tmp);
            urdst[IEN] = urst[IEN] + sqrt_urst_d * sgn_bx * (vbstr - tmp);

        }

        // step 6
        if (s0 >= 0.0) {
            for (int n = 0; n < NVAR; ++n) fw[n] = fl[n];
        } else if (s4 <= 0.0) { 
            for (int n = 0; n < NVAR; ++n) fw[n] = fr[n];
        } else if (s1 >= 0.0) {
            for (int n = 0; n < NVAR; ++n) fw[n] = fl[n] + s0*(ulst[n] - ul[n]);
        } else if (s2 >= 0.0) {
            for (int n = 0; n < NVAR; ++n) fw[n] = fl[n] + s0*(ulst[n] - ul[n]) + s1*(uldst[n] - ulst[n]);
        } else if (s3 > 0.0) {
            for (int n = 0; n < NVAR; ++n) fw[n] = fr[n] + s4*(urst[n] - ur[n]) + s3*(urdst[n] - urst[n]);
        } else {
            for (int n = 0; n < NVAR; ++n) fw[n] = fr[n] + s4*(urst[n] - ur[n]);
        }

        flx[IRHO][k][j][i] = fw[IRHO];
        flx[IVN][k][j][i]  = fw[IV1];
        flx[IVT][k][j][i]  = fw[IV2];
        flx[IVU][k][j][i]  = fw[IV3];
        flx[IEN][k][j][i]  = fw[IEN];
        flx[IBN][k][j][i]  = psim;
        flx[IBT][k][j][i]  = fw[IB2];
        flx[IBU][k][j][i]  = fw[IB3];
        flx[IPSI][k][j][i] = ch*ch*bxm;
    }
}