#include "padd.h"

void PADD_jacob_8x1w_z1eqz2(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    htfe_t t1, t2, t3, t4, t5, t6;

    /* H = x2 - x1. */
    gfp_sub_8x1w(t3, q->x, p->x);

    /* t1 = R = 2 * (y2 - y1). */
    gfp_sub_8x1w(t1, q->y, p->y);
    gfp_add_8x1w(t1, t1, t1);

    /* t2 = HH = H^2. */
    gfp_sqr_8x1w(t2, t3);

    /* t4 = I = 4*HH. */
    gfp_dbl_8x1w(t4, t2);
    gfp_dbl_8x1w(t4, t4);

    /* t5 = J = H * I. */
    gfp_mul_8x1w(t5, t3, t4);

    /* t4 = V = x1 * I. */
    gfp_mul_8x1w(t4, p->x, t4);

    /* x3 = R^2 - J - 2 * V. */
    gfp_sqr_8x1w(r->x, t1);
    gfp_sub_8x1w(r->x, r->x, t5);
    gfp_dbl_8x1w(t6, t4);
    gfp_sub_8x1w(r->x, r->x, t6);

    /* y3 = R * (V - x3) - 2 * Y1 * J. */
    gfp_sub_8x1w(t4, t4, r->x);
    gfp_mul_8x1w(t4, t4, t1);
    gfp_mul_8x1w(t1, p->y, t5);
    gfp_dbl_8x1w(t1, t1);
    gfp_sub_8x1w(r->y, t4, t1);

    /* z3 = 2 * H. */
    gfp_dbl_8x1w(r->z, t3);
}

void PADD_jacob_8x1w_mix(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    htfe_t t0, t1, t2, t3, t4, t5, t6;

    /* t0 = z1^2. */
    gfp_sqr_8x1w(t0, p->z);

    /* t3 = U2 = x2 * z1^2. */
    gfp_mul_8x1w(t3, q->x, t0);

    /* t1 = S2 = y2 * z1^3. */
    gfp_mul_8x1w(t1, t0, p->z);
    gfp_mul_8x1w(t1, t1, q->y);

    /* t3 = H = U2 - x1. */
    gfp_sub_8x1w(t3, t3, p->x);

    /* t1 = R = 2 * (S2 - y1). */
    gfp_sub_8x1w(t1, t1, p->y);
    gfp_dbl_8x1w(t1, t1);

    /* t2 = HH = H^2. */
    gfp_sqr_8x1w(t2, t3);

    /* t4 = I = 4*HH. */
    gfp_dbl_8x1w(t4, t2);
    gfp_dbl_8x1w(t4, t4);

    /* t5 = J = H * I. */
    gfp_mul_8x1w(t5, t3, t4);

    /* t4 = V = x1 * I. */
    gfp_mul_8x1w(t4, p->x, t4);

    /* x3 = R^2 - J - 2 * V. */
    gfp_sqr_8x1w(r->x, t1);
    gfp_sub_8x1w(r->x, r->x, t5);
    gfp_dbl_8x1w(t6, t4);
    gfp_sub_8x1w(r->x, r->x, t6);

    /* y3 = R * (V - x3) - 2 * Y1 * J. */
    gfp_sub_8x1w(t4, t4, r->x);
    gfp_mul_8x1w(t4, t4, t1);
    gfp_mul_8x1w(t1, p->y, t5);
    gfp_dbl_8x1w(t1, t1);
    gfp_sub_8x1w(r->y, t4, t1);

    /* z3 = (z1 + H)^2 - z1^2 - HH. */
    gfp_add_8x1w(r->z, p->z, t3);
    gfp_sqr_8x1w(r->z, r->z);
    gfp_sub_8x1w(r->z, r->z, t0);
    gfp_sub_8x1w(r->z, r->z, t2);
}

void PADD_jacob_8x1w_z1neqz2(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    htfe_t t0, t1, t2, t3, t4, t5, t6;

    /* t0 = z1^2. */
    gfp_sqr_8x1w(t0, p->z);

    /* t1 = z2^2. */
    gfp_sqr_8x1w(t1, q->z);

    /* t2 = U1 = x1 * z2^2. */
    gfp_mul_8x1w(t2, p->x, t1);

    /* t3 = U2 = x2 * z1^2. */
    gfp_mul_8x1w(t3, q->x, t0);

    /* t6 = z1^2 + z2^2. */
    gfp_add_8x1w(t6, t0, t1);

    /* t0 = S2 = y2 * z1^3. */
    gfp_mul_8x1w(t0, t0, p->z);
    gfp_mul_8x1w(t0, t0, q->y);

    /* t1 = S1 = y1 * z2^3. */
    gfp_mul_8x1w(t1, t1, q->z);
    gfp_mul_8x1w(t1, t1, p->y);

    /* t3 = H = U2 - U1. */
    gfp_sub_8x1w(t3, t3, t2);

    /* t0 = R = 2 * (S2 - S1). */
    gfp_sub_8x1w(t0, t0, t1);
    gfp_dbl_8x1w(t0, t0);

    /* t4 = I = (2*H)^2. */
    gfp_dbl_8x1w(t4, t3);
    gfp_sqr_8x1w(t4, t4);

    /* t5 = J = H * I. */
    gfp_mul_8x1w(t5, t3, t4);

    /* t4 = V = U1 * I. */
    gfp_mul_8x1w(t4, t2, t4);

    /* x3 = R^2 - J - 2 * V. */
    gfp_sqr_8x1w(r->x, t0);
    gfp_sub_8x1w(r->x, r->x, t5);
    gfp_dbl_8x1w(t2, t4);
    gfp_sub_8x1w(r->x, r->x, t2);

    /* y3 = R * (V - x3) - 2 * S1 * J. */
    gfp_sub_8x1w(t4, t4, r->x);
    gfp_mul_8x1w(t4, t4, t0);
    gfp_mul_8x1w(t1, t1, t5);
    gfp_dbl_8x1w(t1, t1);
    gfp_sub_8x1w(r->y, t4, t1);

    /* z3 = ((z1 + z2)^2 - z1^2 - z2^2) * H. */
    gfp_add_8x1w(r->z, p->z, q->z);
    gfp_sqr_8x1w(r->z, r->z);
    gfp_sub_8x1w(r->z, r->z, t6);
    gfp_mul_8x1w(r->z, r->z, t3);
}


void PADD_projc_8x1w_z1eqz2(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    htfe_t v, u, uu, vv, R, A;

    /* u = Y2 - Y1. */
    gfp_sub_8x1w(u, q->y, p->y);
    /* v = X2 - X1. */
    gfp_sub_8x1w(v, q->x, p->x);

    /* uu = u ^ 2. */
    gfp_sqr_8x1w(uu, u);
    /* vv = v ^ 2. */
    gfp_sqr_8x1w(vv, v);

    /* Z3 = v ^ 3. */
    gfp_mul_8x1w(r->z, v, vv);
    
    /* R = vv * X1. */
    gfp_mul_8x1w(R, vv, p->x);

    /* A = uu - vvv - 2*R. */
    gfp_sub_8x1w(A, uu, r->z);
    gfp_sub_8x1w(A, A, R);
    gfp_sub_8x1w(A, A, R);

    /* X3 = v * A. */
    gfp_mul_8x1w(r->x, v, A);
    /* Y3 = u * (R - A) - vvv*Y1. */
    gfp_sub_8x1w(R, R, A);
    gfp_mul_8x1w(R, R, u);
    gfp_mul_8x1w(A, p->y, r->z);
    gfp_sub_8x1w(r->y, R, A);

}

void PADD_projc_8x1w_mix(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    // htfe_t t0, t1, t2, t3, t4, t5, t6;

    // gfp_mul_8x1w(t0, p->x, q->x);
    // gfp_mul_8x1w(t1, p->y, q->y);
    // gfp_add_8x1w(t3, q->x, q->y);
    // gfp_add_8x1w(t4, p->x, p->y);
    // gfp_mul_8x1w(t3, t3, t4);
    // gfp_add_8x1w(t4, t0, t1);
    // gfp_sub_8x1w(t3, t3, t4);

    // gfp_mul_8x1w(t4, q->y, p->z);
    // gfp_add_8x1w(t4, t4, p->y);
    // gfp_mul_8x1w(r->y, q->x, p->z);
    // gfp_add_8x1w(r->y, r->y, p->x);
    // /* BLS12-381 b = 3 ep_curve_mul_b3(t2, p->z); */
    // gfp_dbl_8x1w(t2, p->z);
    // gfp_add_8x1w(t2, t2, p->z);
    // gfp_dbl_8x1w(t2, t2);
    // gfp_dbl_8x1w(t2, t2);
    
    // gfp_add_8x1w(r->z, t1, t2);
    // gfp_sub_8x1w(t1, t1, t2);

    // gfp_dbl_8x1w(r->x, t0);
    // gfp_add_8x1w(t0, t0, r->x);
    // /* ep_curve_mul_b3(r->y, r->y);*/
    // gfp_dbl_8x1w(t5, r->y);
    // gfp_add_8x1w(r->y, r->y, t5);
    // gfp_dbl_8x1w(r->y, r->y);
    // gfp_dbl_8x1w(r->y, r->y);

    // gfp_mul_8x1w(r->x, t4, r->y);
    // gfp_mul_8x1w(t2, t3, t1);
    // gfp_sub_8x1w(r->x, t2, r->x);
    // gfp_mul_8x1w(r->y, t0, r->y);
    // gfp_mul_8x1w(t1, t1, r->z);
    // gfp_add_8x1w(r->y, t1, r->y);
    // gfp_mul_8x1w(t0, t0, t3);
    // gfp_mul_8x1w(r->z, r->z, t4);
    // gfp_add_8x1w(r->z, r->z, t0);


    htfe_t v, u, uu, vv, vvv, R, A;

    /* u = Y2*Z1 - Y1. */
    gfp_mul_8x1w(u, q->y, p->z);
    gfp_sub_8x1w(u, u, p->y);
    /* v = X2*Z1 - X1. */
    gfp_mul_8x1w(v, q->x, p->z);
    gfp_sub_8x1w(v, v, p->x);

    /* uu = u ^ 2. */
    gfp_sqr_8x1w(uu, u);
    /* vv = v ^ 2. */
    gfp_sqr_8x1w(vv, v);
    /* vvv = v * vv. */
    gfp_mul_8x1w(vvv, v, vv);

    /* Z3 = vvv * Z1. */
    gfp_mul_8x1w(r->z, vvv, p->z);
    
    /* R = vv * X1. */
    gfp_mul_8x1w(R, vv, p->x);

    /* A = uu * Z1 - vvv - 2*R. */
    gfp_mul_8x1w(A, uu, p->z);
    gfp_sub_8x1w(A, A, vvv);
    gfp_sub_8x1w(A, A, R);
    gfp_sub_8x1w(A, A, R);

    /* X3 = v * A. */
    gfp_mul_8x1w(r->x, v, A);
    /* Y3 = u * (R - A) - vvv*Y1. */
    gfp_sub_8x1w(R, R, A);
    gfp_mul_8x1w(R, R, u);
    gfp_mul_8x1w(A, p->y, vvv);
    gfp_sub_8x1w(r->y, R, A);

}

void PADD_projc_8x1w_z1neqz2(htpoint_t r, const htpoint_t p, const htpoint_t q)
{
    htfe_t t0, t1, t2, t3, t4, t5, t6;

    gfp_mul_8x1w(t0, p->x, q->x);
    gfp_mul_8x1w(t1, p->y, q->y);
    gfp_mul_8x1w(t2, p->z, q->z);
    gfp_add_8x1w(t3, p->x, p->y);
    gfp_add_8x1w(t4, q->x, q->y);
    gfp_mul_8x1w(t3, t3, t4);
    gfp_add_8x1w(t4, t0, t1);
    gfp_sub_8x1w(t3, t3, t4);

    gfp_add_8x1w(t4, p->y, p->z);
    gfp_add_8x1w(t5, q->y, q->z);
    gfp_mul_8x1w(t4, t4, t5);
    gfp_add_8x1w(t5, t1, t2);
    gfp_sub_8x1w(t4, t4, t5);
    gfp_add_8x1w(r->y, q->x, q->z);
    gfp_add_8x1w(r->x, p->x, p->z);
    gfp_mul_8x1w(r->x, r->x, r->y);
    gfp_add_8x1w(r->y, t0, t2);
    gfp_sub_8x1w(r->y, r->x, r->y);
    gfp_dbl_8x1w(r->x, t0);
    gfp_add_8x1w(t0, t0, r->x);

    /* ep_curve_mul_b3(t2, t2);*/
    gfp_dbl_8x1w(t5, t2);
    gfp_add_8x1w(t2, t2, t5);
    gfp_dbl_8x1w(t2, t2);
    gfp_dbl_8x1w(t2, t2);


    gfp_add_8x1w(r->z, t1, t2);
    gfp_sub_8x1w(t1, t1, t2);

    /* ep_curve_mul_b3(r->y, r->y);*/
    gfp_dbl_8x1w(t5, r->y);
    gfp_add_8x1w(r->y, r->y, t5);
    gfp_dbl_8x1w(r->y, r->y);
    gfp_dbl_8x1w(r->y, r->y);

    gfp_mul_8x1w(r->x, t4, r->y);
    gfp_mul_8x1w(t2, t3, t1);
    gfp_sub_8x1w(r->x, t2, r->x);
    gfp_mul_8x1w(r->y, t0, r->y);
    gfp_mul_8x1w(t1, t1, r->z);
    gfp_add_8x1w(r->y, t1, r->y);
    gfp_mul_8x1w(t0, t0, t3);
    gfp_mul_8x1w(r->z, r->z, t4);
    gfp_add_8x1w(r->z, r->z, t0);
}

void PADD_G2_projc_8x1w_z1eqz2(ht2point_t r, const ht2point_t p, const ht2point_t q)
{
    htfe2_t v, u, uu, vv, R, A;
    
    /* u = Y2 - Y1. */
    gfp2_sub_8x1w(u, q->y, p->y);
    /* v = X2 - X1. */
    gfp2_sub_8x1w(v, q->x, p->x);

    /* uu = u ^ 2. */
    gfp2_sqr_8x1w(uu, u);
    /* vv = v ^ 2. */
    gfp2_sqr_8x1w(vv, v);

    /* Z3 = v ^ 3. */
    gfp2_mul_8x1w(r->z, v, vv);
    
    /* R = vv * X1. */
    gfp2_mul_8x1w(R, vv, p->x);

    /* A = uu - vvv - 2*R. */
    gfp2_sub_8x1w(A, uu, r->z);
    gfp2_sub_8x1w(A, A, R);
    gfp2_sub_8x1w(A, A, R);

    /* X3 = v * A. */
    gfp2_mul_8x1w(r->x, v, A);
    /* Y3 = u * (R - A) - vvv*Y1. */
    gfp2_sub_8x1w(R, R, A);
    gfp2_mul_8x1w(R, R, u);
    gfp2_mul_8x1w(A, p->y, r->z);
    gfp2_sub_8x1w(r->y, R, A);
/*----------------------------------------------------*/
    // htfe2_t t0, t1, t2, t3, t4, t5, t6;
    // /* H = x2 - x1. */
    // gfp2_sub_8x1w(t3, q->x, p->x);

    // /* t1 = R = 2 * (y2 - y1). */
    // gfp2_sub_8x1w(t1, q->y, p->y);
    // /* t2 = HH = H^2. */
    // gfp2_sqr_8x1w(t2, t3);
    // /* t5 = J = H * HH. */
    // gfp2_mul_8x1w(t5, t3, t2);

    // /* t4 = V = x1 * HH. */
    // gfp2_mul_8x1w(t4, p->x, t2);

    // /* x3 = R^2 - J - 2 * V. */
    // gfp2_sqr_8x1w(r->x, t1);
    // gfp2_sub_8x1w(r->x, r->x, t5);
    // gfp2_dbl_8x1w(t6, t4);
    // gfp2_sub_8x1w(r->x, r->x, t6);

    // /* y3 = R * (V - x3) - Y1 * J. */
    // gfp2_sub_8x1w(t4, t4, r->x);
    // gfp2_mul_8x1w(t4, t4, t1);
    // gfp2_mul_8x1w(t1, p->y, t5);
    // gfp2_sub_8x1w(r->y, t4, t1);

    // /* z3 = H. */
    // gfp2_copy_8x1w(r->z, t3);
/*-----------------------------------------------------------*/
    // htfe2_t t1, t2, t3, t4, t5, t6;

    // /* H = x2 - x1. */
    // gfp2_sub_8x1w(t3, q->x, p->x);

    // /* t1 = R = 2 * (y2 - y1). */
    // gfp2_sub_8x1w(t1, q->y, p->y);
    // gfp2_add_8x1w(t1, t1, t1);

    // /* t2 = HH = H^2. */
    // gfp2_sqr_8x1w(t2, t3);

    // /* t4 = I = 4*HH. */
    // gfp2_dbl_8x1w(t4, t2);
    // gfp2_dbl_8x1w(t4, t4);

    // /* t5 = J = H * I. */
    // gfp2_mul_8x1w(t5, t3, t4);

    // /* t4 = V = x1 * I. */
    // gfp2_mul_8x1w(t4, p->x, t4);

    // /* x3 = R^2 - J - 2 * V. */
    // gfp2_sqr_8x1w(r->x, t1);
    // gfp2_sub_8x1w(r->x, r->x, t5);
    // gfp2_dbl_8x1w(t6, t4);
    // gfp2_sub_8x1w(r->x, r->x, t6);

    // /* y3 = R * (V - x3) - 2 * Y1 * J. */
    // gfp2_sub_8x1w(t4, t4, r->x);
    // gfp2_mul_8x1w(t4, t4, t1);
    // gfp2_mul_8x1w(t1, p->y, t5);
    // gfp2_dbl_8x1w(t1, t1);
    // gfp2_sub_8x1w(r->y, t4, t1);

    // /* z3 = 2 * H. */
    // gfp2_dbl_8x1w(r->z, t3);

}
void PADD_G2_projc_8x1w_mix(ht2point_t r, const ht2point_t p, const ht2point_t q)
{
    // htfe2_t t0, t1, t2, t3, t4, t5, t6;

    // gfp2_mul_8x1w(t0, p->x, q->x);
    // gfp2_mul_8x1w(t1, p->y, q->y);
    // gfp2_add_8x1w(t3, q->x, q->y);
    // gfp2_add_8x1w(t4, p->x, p->y);
    // gfp2_mul_8x1w(t3, t3, t4);
    // gfp2_add_8x1w(t4, t0, t1);
    // gfp2_sub_8x1w(t3, t3, t4);

    // gfp2_mul_8x1w(t4, q->y, p->z);
    // gfp2_add_8x1w(t4, t4, p->y);
    // gfp2_mul_8x1w(r->y, q->x, p->z);
    // gfp2_add_8x1w(r->y, r->y, p->x);
    // /* BLS12-381 b = 3 ep_curve_mul_b3(t2, p->z); */
    // gfp2_mul_dig_8x1w(p->z, p->z);
    // gfp2_dbl_8x1w(t2, p->z);
    // gfp2_add_8x1w(t2, t2, p->z);
    // gfp2_dbl_8x1w(t2, t2);
    // gfp2_dbl_8x1w(t2, t2);
    
    // gfp2_add_8x1w(r->z, t1, t2);
    // gfp2_sub_8x1w(t1, t1, t2);

    // gfp2_dbl_8x1w(r->x, t0);
    // gfp2_add_8x1w(t0, t0, r->x);
    // /* ep_curve_mul_b3(r->y, r->y);*/
    // gfp2_mul_dig_8x1w(r->y, r->y);
    // gfp2_dbl_8x1w(t5, r->y);
    // gfp2_add_8x1w(r->y, r->y, t5);
    // gfp2_dbl_8x1w(r->y, r->y);
    // gfp2_dbl_8x1w(r->y, r->y);

    // gfp2_mul_8x1w(r->x, t4, r->y);
    // gfp2_mul_8x1w(t2, t3, t1);
    // gfp2_sub_8x1w(r->x, t2, r->x);
    // gfp2_mul_8x1w(r->y, t0, r->y);
    // gfp2_mul_8x1w(t1, t1, r->z);
    // gfp2_add_8x1w(r->y, t1, r->y);
    // gfp2_mul_8x1w(t0, t0, t3);
    // gfp2_mul_8x1w(r->z, r->z, t4);
    // gfp2_add_8x1w(r->z, r->z, t0);    

    /*------------------------------------------------------------*/
    htfe2_t v, u, uu, vv, vvv, R, A;

    /* u = Y2*Z1 - Y1. */
    gfp2_mul_8x1w(u, q->y, p->z);
    gfp2_sub_8x1w(u, u, p->y);
    /* v = X2*Z1 - X1. */
    gfp2_mul_8x1w(v, q->x, p->z);
    gfp2_sub_8x1w(v, v, p->x);

    /* uu = u ^ 2. */
    gfp2_sqr_8x1w(uu, u);
    /* vv = v ^ 2. */
    gfp2_sqr_8x1w(vv, v);
    /* vvv = v * vv. */
    gfp2_mul_8x1w(vvv, v, vv);

    /* Z3 = vvv * Z1. */
    gfp2_mul_8x1w(r->z, vvv, p->z);
    
    /* R = vv * X1. */
    gfp2_mul_8x1w(R, vv, p->x);

    /* A = uu * Z1 - vvv - 2*R. */
    gfp2_mul_8x1w(A, uu, p->z);
    gfp2_sub_8x1w(A, A, vvv);
    gfp2_sub_8x1w(A, A, R);
    gfp2_sub_8x1w(A, A, R);

    /* X3 = v * A. */
    gfp2_mul_8x1w(r->x, v, A);
    /* Y3 = u * (R - A) - vvv*Y1. */
    gfp2_sub_8x1w(R, R, A);
    gfp2_mul_8x1w(R, R, u);
    gfp2_mul_8x1w(A, p->y, vvv);
    gfp2_sub_8x1w(r->y, R, A);
}
void PADD_G2_projc_8x1w_z1neqz2(ht2point_t r, const ht2point_t p, const ht2point_t q)
{
    htfe2_t t0, t1, t2, t3, t4, t5, t6;

    gfp2_mul_8x1w(t0, p->x, q->x);
    gfp2_mul_8x1w(t1, p->y, q->y);
    gfp2_mul_8x1w(t2, p->z, q->z);
    gfp2_add_8x1w(t3, p->x, p->y);
    gfp2_add_8x1w(t4, q->x, q->y);
    gfp2_mul_8x1w(t3, t3, t4);
    gfp2_add_8x1w(t4, t0, t1);
    gfp2_sub_8x1w(t3, t3, t4);

    gfp2_add_8x1w(t4, p->y, p->z);
    gfp2_add_8x1w(t5, q->y, q->z);
    gfp2_mul_8x1w(t4, t4, t5);
    gfp2_add_8x1w(t5, t1, t2);
    gfp2_sub_8x1w(t4, t4, t5);
    gfp2_add_8x1w(r->y, q->x, q->z);
    gfp2_add_8x1w(r->x, p->x, p->z);
    gfp2_mul_8x1w(r->x, r->x, r->y);
    gfp2_add_8x1w(r->y, t0, t2);
    gfp2_sub_8x1w(r->y, r->x, r->y);
    gfp2_dbl_8x1w(r->x, t0);
    gfp2_add_8x1w(t0, t0, r->x);

    /* ep_curve_mul_b3(t2, t2);*/
    gfp2_mul_dig_8x1w(t2, t2);
    gfp2_dbl_8x1w(t5, t2);
    gfp2_add_8x1w(t2, t2, t5);
    gfp2_dbl_8x1w(t2, t2);
    gfp2_dbl_8x1w(t2, t2);


    gfp2_add_8x1w(r->z, t1, t2);
    gfp2_sub_8x1w(t1, t1, t2);

    /* ep_curve_mul_b3(r->y, r->y);*/
    gfp2_mul_dig_8x1w(r->y, r->y);
    gfp2_dbl_8x1w(t5, r->y);
    gfp2_add_8x1w(r->y, r->y, t5);
    gfp2_dbl_8x1w(r->y, r->y);
    gfp2_dbl_8x1w(r->y, r->y);

    gfp2_mul_8x1w(r->x, t4, r->y);
    gfp2_mul_8x1w(t2, t3, t1);
    gfp2_sub_8x1w(r->x, t2, r->x);
    gfp2_mul_8x1w(r->y, t0, r->y);
    gfp2_mul_8x1w(t1, t1, r->z);
    gfp2_add_8x1w(r->y, t1, r->y);
    gfp2_mul_8x1w(t0, t0, t3);
    gfp2_mul_8x1w(r->z, r->z, t4);
    gfp2_add_8x1w(r->z, r->z, t0);    
}