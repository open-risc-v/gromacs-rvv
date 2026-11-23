/*
 * RVV/vector kernel for the GROMACS SETTLE constraint algorithm.
 *
 * The scalar reference lives in mdlib/settle.cpp (settleTemplate, packSize==1).
 * SETTLE is per-molecule: no accumulation crosses molecules except the virial
 * sum, so binding one lane to one molecule keeps every floating-point
 * operation sequence elementwise identical to the scalar code.
 *
 * Lanes = 8 floats (32 B vector).  Coordinates are AoS (x[i*3+d]), so a lane
 * block is gathered with stride 3 and scattered back the same way.
 *
 * Functions that are not plain +-* (gmx::invsqrt, gmx::fma, max) are applied
 * lane by lane through the *scalar* implementations, so their rounding
 * behaviour is exactly the reference one.  The virial sum is accumulated
 * in molecule order (lane 0..n) rather than reduced across lanes, again to
 * keep the addition order identical.
 *
 * Only instantiated for the packSize==1 path (GMX_SIMD=None).  PBC folding
 * is delegated to the scalar pbc_dx_aiuc, one call per lane.
 */
#ifndef GROMACS_SETTLE_RVV_H
#define GROMACS_SETTLE_RVV_H

typedef float settle_v8 __attribute__((vector_size(32)));
#define SETTLE_RVV_LANES 8

static inline settle_v8 settle_rvv_bcast(real s)
{
    settle_v8 r;
    for (int i = 0; i < SETTLE_RVV_LANES; i++)
    {
        r[i] = s;
    }
    return r;
}

static inline settle_v8 settle_rvv_invsqrt(settle_v8 a)
{
    settle_v8 r;
    for (int i = 0; i < SETTLE_RVV_LANES; i++)
    {
        r[i] = gmx::invsqrt(a[i]);
    }
    return r;
}

static inline settle_v8 settle_rvv_fma(settle_v8 a, settle_v8 b, settle_v8 c)
{
    settle_v8 r;
    for (int i = 0; i < SETTLE_RVV_LANES; i++)
    {
        r[i] = gmx::fma(a[i], b[i], c[i]);
    }
    return r;
}

static inline settle_v8 settle_rvv_max(settle_v8 a, settle_v8 b)
{
    settle_v8 r;
    for (int i = 0; i < SETTLE_RVV_LANES; i++)
    {
        r[i] = (a[i] > b[i]) ? a[i] : b[i];
    }
    return r;
}

static void settle_rvv(const SettleData& settled,
                       const t_pbc*     pbc,
                       int               settleStart,
                       int               settleEnd,
                       const real*       x,
                       real*             xprime,
                       real              invdt,
                       real*             v,
                       bool              bCalcVirial,
                       tensor            vir_r_m_dr,
                       bool*             bErrorHasOccurred)
{
    const SettleParameters* p = &settled.parametersMassWeighted();

    settle_v8 wh   = settle_rvv_bcast(p->wh);
    settle_v8 rc   = settle_rvv_bcast(p->rc);
    settle_v8 ra   = settle_rvv_bcast(p->ra);
    settle_v8 rb   = settle_rvv_bcast(p->rb);
    settle_v8 irc2 = settle_rvv_bcast(p->irc2);
    settle_v8 mO   = settle_rvv_bcast(p->mO);
    settle_v8 mH   = settle_rvv_bcast(p->mH);
    settle_v8 almost_zero = settle_rvv_bcast(real(1e-12));

    bool bError = false;
    real sum_r_m_dr[DIM][DIM];
    for (int d2 = 0; d2 < DIM; d2++)
    {
        for (int d = 0; d < DIM; d++)
        {
            sum_r_m_dr[d2][d] = 0;
        }
    }

    for (int i = settleStart; i < settleEnd; i += SETTLE_RVV_LANES)
    {
        int nl = settleEnd - i;
        if (nl > SETTLE_RVV_LANES)
        {
            nl = SETTLE_RVV_LANES;
        }

        settle_v8 xo[DIM], xh2[DIM], xh3[DIM];
        settle_v8 xpo[DIM], xph2[DIM], xph3[DIM];
        settle_v8 vo[DIM], vh2[DIM], vh3[DIM];
        for (int d = 0; d < DIM; d++)
        {
            xo[d] = xh2[d] = xh3[d] = xpo[d] = xph2[d] = xph3[d] = settle_rvv_bcast(0);
            vo[d] = vh2[d] = vh3[d] = settle_rvv_bcast(0);
        }
        for (int l = 0; l < nl; l++)
        {
            const int io  = settled.ow1()[i + l] * DIM;
            const int ih2 = settled.hw2()[i + l] * DIM;
            const int ih3 = settled.hw3()[i + l] * DIM;
            for (int d = 0; d < DIM; d++)
            {
                xo[d][l]   = x[io + d];
                xh2[d][l]  = x[ih2 + d];
                xh3[d][l]  = x[ih3 + d];
                xpo[d][l]  = xprime[io + d];
                xph2[d][l] = xprime[ih2 + d];
                xph3[d][l] = xprime[ih3 + d];
            }
            if (v != nullptr)
            {
                for (int d = 0; d < DIM; d++)
                {
                    vo[d][l]  = v[io + d];
                    vh2[d][l] = v[ih2 + d];
                    vh3[d][l] = v[ih3 + d];
                }
            }
        }

        /* Distances go through the *scalar* pbc_dx_aiuc per lane: the PBC
         * folding must match the reference bit for bit, so it is not
         * reimplemented with masks here.
         */
        settle_v8 dist21[DIM], dist31[DIM], doh2[DIM], doh3[DIM];
        for (int d = 0; d < DIM; d++)
        {
            dist21[d] = dist31[d] = doh2[d] = doh3[d] = settle_rvv_bcast(0);
        }
        for (int l = 0; l < nl; l++)
        {
            rvec a, b, dx;
            for (int d = 0; d < DIM; d++)
            {
                a[d] = xh2[d][l];
                b[d] = xo[d][l];
            }
            pbc_dx_aiuc(pbc, a, b, dx);
            for (int d = 0; d < DIM; d++)
            {
                dist21[d][l] = dx[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                a[d] = xh3[d][l];
            }
            pbc_dx_aiuc(pbc, a, b, dx);
            for (int d = 0; d < DIM; d++)
            {
                dist31[d][l] = dx[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                a[d] = xph2[d][l];
                b[d] = xpo[d][l];
            }
            pbc_dx_aiuc(pbc, a, b, dx);
            for (int d = 0; d < DIM; d++)
            {
                doh2[d][l] = dx[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                a[d] = xph3[d][l];
            }
            pbc_dx_aiuc(pbc, a, b, dx);
            for (int d = 0; d < DIM; d++)
            {
                doh3[d][l] = dx[d];
            }
        }

        settle_v8 a1[DIM], b1[DIM], c1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            a1[d] = -(doh2[d] + doh3[d]) * wh;
            b1[d] = doh2[d] + a1[d];
            c1[d] = doh3[d] + a1[d];
        }

        settle_v8 xakszd = dist21[YY] * dist31[ZZ] - dist21[ZZ] * dist31[YY];
        settle_v8 yakszd = dist21[ZZ] * dist31[XX] - dist21[XX] * dist31[ZZ];
        settle_v8 zakszd = dist21[XX] * dist31[YY] - dist21[YY] * dist31[XX];
        settle_v8 xaksxd = a1[YY] * zakszd - a1[ZZ] * yakszd;
        settle_v8 yaksxd = a1[ZZ] * xakszd - a1[XX] * zakszd;
        settle_v8 zaksxd = a1[XX] * yakszd - a1[YY] * xakszd;
        settle_v8 xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        settle_v8 yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        settle_v8 zaksyd = xakszd * yaksxd - yakszd * xaksxd;

        settle_v8 axlng = settle_rvv_invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        settle_v8 aylng = settle_rvv_invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        settle_v8 azlng = settle_rvv_invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        settle_v8 trns1[DIM], trns2[DIM], trns3[DIM];
        trns1[XX] = xaksxd * axlng;
        trns2[XX] = yaksxd * axlng;
        trns3[XX] = zaksxd * axlng;
        trns1[YY] = xaksyd * aylng;
        trns2[YY] = yaksyd * aylng;
        trns3[YY] = zaksyd * aylng;
        trns1[ZZ] = xakszd * azlng;
        trns2[ZZ] = yakszd * azlng;
        trns3[ZZ] = zakszd * azlng;

        settle_v8 b0d[2], c0d[2];
        for (int d = 0; d < 2; d++)
        {
            b0d[d] = trns1[d] * dist21[XX] + trns2[d] * dist21[YY] + trns3[d] * dist21[ZZ];
            c0d[d] = trns1[d] * dist31[XX] + trns2[d] * dist31[YY] + trns3[d] * dist31[ZZ];
        }

        settle_v8 a1d_z = trns1[ZZ] * a1[XX] + trns2[ZZ] * a1[YY] + trns3[ZZ] * a1[ZZ];
        settle_v8 b1d[DIM], c1d[DIM];
        for (int d = 0; d < DIM; d++)
        {
            b1d[d] = trns1[d] * b1[XX] + trns2[d] * b1[YY] + trns3[d] * b1[ZZ];
            c1d[d] = trns1[d] * c1[XX] + trns2[d] * c1[YY] + trns3[d] * c1[ZZ];
        }

        settle_v8 sinphi = a1d_z * settle_rvv_invsqrt(ra * ra);
        settle_v8 tmp2   = settle_rvv_bcast(real(1)) - sinphi * sinphi;
        for (int l = 0; l < nl; l++)
        {
            if (tmp2[l] <= real(1e-12))
            {
                bError = true;
            }
        }
        tmp2 = settle_rvv_max(tmp2, almost_zero);
        settle_v8 tmp    = settle_rvv_invsqrt(tmp2);
        settle_v8 cosphi = tmp2 * tmp;
        settle_v8 sinpsi = (b1d[ZZ] - c1d[ZZ]) * irc2 * tmp;
        tmp2             = settle_rvv_bcast(real(1)) - sinpsi * sinpsi;
        settle_v8 cospsi = tmp2 * settle_rvv_invsqrt(tmp2);

        settle_v8 a2d_y = ra * cosphi;
        settle_v8 b2d_x = -rc * cospsi;
        settle_v8 t1    = -rb * cosphi;
        settle_v8 t2    = rc * sinpsi * sinphi;
        settle_v8 b2d_y = t1 - t2;
        settle_v8 c2d_y = t1 + t2;

        settle_v8 alpha = b2d_x * (b0d[XX] - c0d[XX]) + b0d[YY] * b2d_y + c0d[YY] * c2d_y;
        settle_v8 beta  = b2d_x * (c0d[YY] - b0d[YY]) + b0d[XX] * b2d_y + c0d[XX] * c2d_y;
        settle_v8 gamma = b0d[XX] * b1d[YY] - b1d[XX] * b0d[YY] + c0d[XX] * c1d[YY] - c1d[XX] * c0d[YY];
        settle_v8 al2be2 = alpha * alpha + beta * beta;
        tmp2             = al2be2 - gamma * gamma;
        settle_v8 sinthe = (alpha * gamma - beta * tmp2 * settle_rvv_invsqrt(tmp2))
                           * settle_rvv_invsqrt(al2be2 * al2be2);

        tmp2             = settle_rvv_bcast(real(1)) - sinthe * sinthe;
        settle_v8 costhe = tmp2 * settle_rvv_invsqrt(tmp2);

        settle_v8 a3d[DIM], b3d[DIM], c3d[DIM];
        a3d[XX] = -a2d_y * sinthe;
        a3d[YY] = a2d_y * costhe;
        a3d[ZZ] = a1d_z;
        b3d[XX] = b2d_x * costhe - b2d_y * sinthe;
        b3d[YY] = b2d_x * sinthe + b2d_y * costhe;
        b3d[ZZ] = b1d[ZZ];
        c3d[XX] = -b2d_x * costhe - c2d_y * sinthe;
        c3d[YY] = -b2d_x * sinthe + c2d_y * costhe;
        c3d[ZZ] = c1d[ZZ];

        settle_v8 a3[DIM], b3[DIM], c3[DIM];
        a3[XX] = trns1[XX] * a3d[XX] + trns1[YY] * a3d[YY] + trns1[ZZ] * a3d[ZZ];
        a3[YY] = trns2[XX] * a3d[XX] + trns2[YY] * a3d[YY] + trns2[ZZ] * a3d[ZZ];
        a3[ZZ] = trns3[XX] * a3d[XX] + trns3[YY] * a3d[YY] + trns3[ZZ] * a3d[ZZ];
        b3[XX] = trns1[XX] * b3d[XX] + trns1[YY] * b3d[YY] + trns1[ZZ] * b3d[ZZ];
        b3[YY] = trns2[XX] * b3d[XX] + trns2[YY] * b3d[YY] + trns2[ZZ] * b3d[ZZ];
        b3[ZZ] = trns3[XX] * b3d[XX] + trns3[YY] * b3d[YY] + trns3[ZZ] * b3d[ZZ];
        c3[XX] = trns1[XX] * c3d[XX] + trns1[YY] * c3d[YY] + trns1[ZZ] * c3d[ZZ];
        c3[YY] = trns2[XX] * c3d[XX] + trns2[YY] * c3d[YY] + trns2[ZZ] * c3d[ZZ];
        c3[ZZ] = trns3[XX] * c3d[XX] + trns3[YY] * c3d[YY] + trns3[ZZ] * c3d[ZZ];

        settle_v8 dxOw1[DIM], dxHw2[DIM], dxHw3[DIM];
        for (int d = 0; d < DIM; d++)
        {
            dxOw1[d] = a3[d] - a1[d];
            dxHw2[d] = b3[d] - b1[d];
            dxHw3[d] = c3[d] - c1[d];
        }
        for (int d = 0; d < DIM; d++)
        {
            xpo[d]  = xpo[d] + dxOw1[d];
            xph2[d] = xph2[d] + dxHw2[d];
            xph3[d] = xph3[d] + dxHw3[d];
        }

        for (int l = 0; l < nl; l++)
        {
            const int io  = settled.ow1()[i + l] * DIM;
            const int ih2 = settled.hw2()[i + l] * DIM;
            const int ih3 = settled.hw3()[i + l] * DIM;
            for (int d = 0; d < DIM; d++)
            {
                xprime[io + d]  = xpo[d][l];
                xprime[ih2 + d] = xph2[d][l];
                xprime[ih3 + d] = xph3[d][l];
            }
        }

        if (v != nullptr)
        {
            for (int d = 0; d < DIM; d++)
            {
                vo[d]  = settle_rvv_fma(dxOw1[d], settle_rvv_bcast(invdt), vo[d]);
                vh2[d] = settle_rvv_fma(dxHw2[d], settle_rvv_bcast(invdt), vh2[d]);
                vh3[d] = settle_rvv_fma(dxHw3[d], settle_rvv_bcast(invdt), vh3[d]);
            }
            for (int l = 0; l < nl; l++)
            {
                const int io  = settled.ow1()[i + l] * DIM;
                const int ih2 = settled.hw2()[i + l] * DIM;
                const int ih3 = settled.hw3()[i + l] * DIM;
                for (int d = 0; d < DIM; d++)
                {
                    v[io + d]  = vo[d][l];
                    v[ih2 + d] = vh2[d][l];
                    v[ih3 + d] = vh3[d][l];
                }
            }
        }

        if (bCalcVirial)
        {
            settle_v8 filter = settle_rvv_bcast(0);
            for (int l = 0; l < nl; l++)
            {
                filter[l] = settled.virfac()[i + l];
            }
            settle_v8 mOf = filter * mO;
            settle_v8 mHf = filter * mH;
            settle_v8 mdo[DIM], mdb[DIM], mdc[DIM];
            for (int d = 0; d < DIM; d++)
            {
                mdb[d] = mHf * dxHw2[d];
                mdc[d] = mHf * dxHw3[d];
                mdo[d] = mOf * dxOw1[d] + mdb[d] + mdc[d];
            }
            for (int l = 0; l < nl; l++)
            {
                for (int d2 = 0; d2 < DIM; d2++)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        sum_r_m_dr[d2][d] = sum_r_m_dr[d2][d]
                                - (xo[d2][l] * mdo[d][l] + dist21[d2][l] * mdb[d][l]
                                   + dist31[d2][l] * mdc[d][l]);
                    }
                }
            }
        }
    }

    if (bCalcVirial)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            for (int d = 0; d < DIM; d++)
            {
                vir_r_m_dr[d2][d] += sum_r_m_dr[d2][d];
            }
        }
    }

    if (bError)
    {
        *bErrorHasOccurred = true;
    }
}

#endif /* GROMACS_SETTLE_RVV_H */
