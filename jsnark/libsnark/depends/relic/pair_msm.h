#ifndef PAIR_MSM_H
#define PAIR_MSM_H

#include "include/relic.h"
#include "data_struct.h"
#include "types.h"

typedef struct
{
    int NWINS;
    int WBITS;
    int BUCKETNUM;
    pair_Bucket* bucket;

    pair_buf1_st buf1[AVX_WAY];
    int cnt1;
    pair_buf2_st buf2[AVX_WAY];
    int cnt2;

    pair52 *avx3_out[AVX_WAY];
    int avx3_cnt;
    
    pair_circle_queue T;
;
} pair_msm_ctx;


// void Pippenger_pair_old(int m, pair G, const bn_t *k, const pair *P);
void Pippenger_pair_old(int m, pair G, const bn_t *k, const pair *P, pair_msm_ctx *ctx);
void MSM_pair_old(int m, pair G, const bn_t *k, const pair *P);
void init_pair_ctx(pair_msm_ctx *ctx);
void avx2_pair(pair_msm_ctx *ctx, int eff_times);
void avx1_pair(pair_msm_ctx *ctx, int eff_times);
void put_pair_into_bucket(const pair *P, const uint32_t b_id, int p_id, pair_msm_ctx *ctx);
void pippenger_pair_fast(const bn_t *k, const pair *P, int num, pair result, pair *M, pair_msm_ctx *ctx);


#endif