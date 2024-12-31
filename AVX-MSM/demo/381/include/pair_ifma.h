#ifndef PAIR_IFMA_H
#define PAIR_IFMA_H

#include "pair_msm.h"
#include "padd.h"
#include "utils.h"

/** PADD options. */
#define PROJC_eq     4
#define PROJC_mix    5
#define PROJC_neq    6

void ep2_norm_projc(ep2_t r, const ep2_t p);
void fp2_mul_1i(fp2_t r, fp2_t  a);
void ep2_add_projective_imp(ep2_t r, const ep2_t p, const ep2_t q);
void ep2_add_projective(ep2_t r, const ep2_t p, const ep2_t q);
void ep2_dbl_projective_imp(ep2_t r, const ep2_t p);
void ep2_dbl_projective(ep2_t r, const ep2_t p);


void simdPADDaff_G1(uint64_t u1x[][HT_NWORDS], uint64_t u1y[][HT_NWORDS], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], pair52 p3[8]);
void simdPADDaff_G2(uint64_t u1x[][2][HT_NWORDS], uint64_t u1y[][2][HT_NWORDS], uint64_t u2x[][2][HT_NWORDS], uint64_t u2y[][2][HT_NWORDS], pair52 p3[8]);
void simdPADDprj_G1(pair_buf2_st* buf, pair52 p3[8]);
void simdPADDprj_G2(pair_buf2_st* buf, pair52 p3[8]);
void simdPADDmix_G1(pair52_st *bucket[AVX_WAY], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], int num);
void simdPADDmix_G2(pair52_st *bucket[AVX_WAY], uint64_t u2x[][2][HT_NWORDS], uint64_t u2y[][2][HT_NWORDS], int num);
buf2_st* check_eq_pair_id(pair_msm_ctx *ctx, int id);
void avx2_pair_withIFMA(pair_msm_ctx *ctx, int eff_times);
void pair_state_trans_T_withIFMA(pair_msm_ctx *ctx);
void avx1_pair_withIFMA(pair_msm_ctx *ctx, int eff_times);
void put_pair_into_bucket_withIFMA(const pair *P, const uint32_t b_id, int p_id, pair_msm_ctx *ctx);
void pippenger_pair_ifma(const bn_t *k, const pair *P, int num, pair result, pair *M, pair_msm_ctx *ctx);
void pair_mont64_to_mont52(pair52 q, pair p);
void pair_mont52_to_mont64(pair52 q, pair p);


#endif // PAIR_IFMA_H