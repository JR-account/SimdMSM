#ifndef _GFP2ARIT_H
#define _GFP2ARIT_H

#include "gfparith.h"

typedef htfe_t htfe2_t[2];

void gfp2_add_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b);
void gfp2_dbl_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_sub_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b);
void gfp2_mul_8x1w(htfe2_t r, const htfe2_t a, const htfe2_t b);
void gfp2_sqr_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_rdcp_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_mul_dig_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_zero_8x1w(htfe2_t r);
void gfp2_mont_relic2avx_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_mont_avx2relic_8x1w(htfe2_t r, const htfe2_t a);
void gfp2_copy_8x1w(htfe2_t r, const htfe2_t a);

// void gfp_num2mont_8x1w(htfe_t r, const htfe_t a);
// void gfp_mont2num_8x1w(htfe_t r, const htfe_t a);
// int  gfp_iszero_8x1w(const htfe_t a);

#endif
