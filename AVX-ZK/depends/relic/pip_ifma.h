#ifndef PIP_IFMA_H
#define PIP_IFMA_H

#include "pippenger.h"
#include "padd.h"
#include "utils.h"

/** PADD options. */
#define JABOC_eq     1
#define JABOC_mix    2
#define JABOC_neq    3
#define PROJC_eq     4
#define PROJC_mix    5
#define PROJC_neq    6

void simdPADDaff(uint64_t u1x[][HT_NWORDS], uint64_t u1y[][HT_NWORDS], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], point52_t p3[8]);
void simdPADDprj(buf2_st* buf, point52_t p3[8]);
void simdPADDmix(point52 *bucket[AVX_WAY], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], int num);

int point52_is_infty(const point52_t p);
void avx2_withIFMA(pip_ctx *ctx, int eff_times);
void state_trans_T_withIFMA(pip_ctx *ctx);
void avx1_withIFMA(pip_ctx *ctx, int eff_times);
void put_into_bucket_withIFMA(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx);
void pippenger_ifma(const bn_t *k, const ep_t *P, int num, ep_t result, ep_t *M, pip_ctx *ctx);
void mont64_to_mont52(point52_t q, ep_t p);
void mont52_to_mont64(point52_t q, ep_t p);
#endif // PIP_IFMA_H