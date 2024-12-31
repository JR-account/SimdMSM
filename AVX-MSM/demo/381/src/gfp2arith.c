#include "gfp2arith.h"

void gfp2_add_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b)
{
    gfp_add_8x1w(r[0], a[0], b[0]);
    gfp_add_8x1w(r[1], a[1], b[1]);
}

void gfp2_dbl_8x1w(htfe2_t r, const htfe2_t a)
{
    gfp_dbl_8x1w(r[0], a[0]);
    gfp_dbl_8x1w(r[1], a[1]);
}

void gfp2_sub_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b)
{
    gfp_sub_8x1w(r[0], a[0], b[0]);
    gfp_sub_8x1w(r[1], a[1], b[1]);
}

void gfp2_mul_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b)
{
    htfe_t t0, t1, t2;
    htfe2_t c;
    gfp_mul_8x1w(t0, a[0], b[0]); // t0 = a0*b0
    gfp_mul_8x1w(t1, a[1], b[1]); // t1 = a1*b1
    gfp_sub_8x1w(c[0], t0, t1);   // r0 = t0 - t1 = a0*b0 - a1*b1
    gfp_add_8x1w(t2, a[0], a[1]); // t2 = a0 + a1
    gfp_add_8x1w(c[1], b[0], b[1]); // r1 = b0 + b1
    gfp_mul_8x1w(c[1], t2, c[1]);   // r1 = t2*r1 = (a0 + a1)*(b0 + b1)
    gfp_sub_8x1w(c[1], c[1], t0);   // r1 = r1 - t0 = (a0 + a1)*(b0 + b1) - a0*b0
    gfp_sub_8x1w(c[1], c[1], t1);   // r1 = r1 - t1 = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1

    gfp_copy_8x1w(r[0], c[0]);
    gfp_copy_8x1w(r[1], c[1]);
}
void gfp2_sqr_8x1w(htfe2_t r, const htfe2_t a)
{
    htfe_t t0;
    htfe2_t c;
    gfp_add_8x1w(t0, a[0], a[1]); // t0 = a0 + a1
    gfp_sub_8x1w(c[0], a[0], a[1]); // r0 = a0 - a1
    gfp_mul_8x1w(c[0], t0, c[0]);   // r0 = t0*r0 = (a0 + a1)*(a0 - a1)
    gfp_mul_8x1w(c[1], a[0], a[1]); // r1 = a0*a1
    gfp_add_8x1w(c[1], c[1], c[1]); // r1 = r1 + r1 = 2*a0*a1

    gfp_copy_8x1w(r[0], c[0]);
    gfp_copy_8x1w(r[1], c[1]);
}

void gfp2_mul_dig_8x1w(htfe2_t r, const htfe2_t a)
{
    htfe2_t c;
    gfp_sub_8x1w(c[0], a[0], a[1]);
    gfp_add_8x1w(c[1], a[0], a[1]);

    gfp_copy_8x1w(r[0], c[0]);
    gfp_copy_8x1w(r[1], c[1]);
}

void gfp2_rdcp_8x1w(htfe2_t r, const htfe2_t a)
{
    gfp_rdcp_8x1w(r[0], a[0]);
    gfp_rdcp_8x1w(r[1], a[1]);
}

void gfp2_zero_8x1w(htfe2_t r)
{
    gfp_zero_8x1w(r[0]);
    gfp_zero_8x1w(r[1]);
}

void gfp2_mont_relic2avx_8x1w(htfe2_t r, const htfe2_t a)
{
    gfp_mont_relic2avx_8x1w(r[0], a[0]);
    gfp_mont_relic2avx_8x1w(r[1], a[1]);
}

void gfp2_mont_avx2relic_8x1w(htfe2_t r, const htfe2_t a)
{
    gfp_mont_avx2relic_8x1w(r[0], a[0]);
    gfp_mont_avx2relic_8x1w(r[1], a[1]);
}

void gfp2_copy_8x1w(htfe2_t r, const htfe2_t a)
{
    gfp_copy_8x1w(r[0], a[0]);
    gfp_copy_8x1w(r[1], a[1]);
}