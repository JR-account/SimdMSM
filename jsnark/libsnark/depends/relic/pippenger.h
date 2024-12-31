#ifndef PIPPENGER_H
#define PIPPENGER_H

#include "types.h"
#include "data_struct.h"
typedef struct
{
    int NWINS;
    int WBITS;
    int BUCKETNUM;

    Bucket* bucket;
    buf1_st buf1[AVX_WAY];
    int cnt1;
    buf2_st buf2[AVX_WAY];
    int cnt2;

    point52_t *avx3_out[AVX_WAY];
    int avx3_cnt;
    circle_queue T;
} pip_ctx;

// uint32_t get_wval(const bn_t d, uint32_t offset);
uint32_t get_wval(const bn_t d, uint32_t offset, int WBITS);

// void Pippenger_old(int m, ep_t G, const bn_t *k, const ep_t *P);
void Pippenger_old(int m, ep_t G, const bn_t *k, const ep_t *P,  pip_ctx *ctx);
void MSM_old(int m, ep_t G, const bn_t *k, const ep_t *P);
void init(pip_ctx *ctx);
void avx2(pip_ctx *ctx, int eff_times);
buf2_st* check_eq_id(pip_ctx *ctx, int id);
void avx1(pip_ctx *ctx, int eff_times);
void put_into_bucket(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx);
void pippenger_fast(const bn_t *k, const ep_t *P, int num, ep_t result, ep_t *M, pip_ctx *ctx);

#endif /* !PIPPENGER_H */