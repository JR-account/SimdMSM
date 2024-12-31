#include "pair_threads.h"
#include <math.h>
#include <omp.h>
#include <time.h>

void pair_window(const bn_t *k, const pair *P, int num, pair Gu, pair_msm_ctx *ctx, int index)
{
    int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

    uint32_t b;
	int i, a;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    uint64_t v2x[AVX_WAY][2][HT_NWORDS], v2y[AVX_WAY][2][HT_NWORDS];
    pair tmp;

    pair *M = (pair*)malloc(ctx->BUCKETNUM * sizeof(pair));

    init_pair_ctx(ctx);
    for (int j = 0; j < num; j++) {
        b = get_wval(k[j], index, WBITS);
        if (b != 0) {
            // ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
            put_pair_into_bucket_withIFMA(P, (uint32_t)b - 1, j, ctx);
        }
    }

    /* deal with the remains in buf*/
    /* remain buf1*/
    if(ctx->cnt1 != 0)
    {
        avx1_pair_withIFMA(ctx, ctx->cnt1);
    }
    /* remain buf2*/
    while((!pair_queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
    {
        avx2_pair_withIFMA(ctx, ctx->cnt2);
        pair_state_trans_T_withIFMA(ctx);
    }

    /* remain buf0*/
    for(a = 0; a < BUCKETNUM; a++)
    {
        if (ctx->bucket[a].wait != NULL)
        {
            // ep_add_projc(ctx->bucket[ctx->dblpoint[a].id].B, ctx->bucket[ctx->dblpoint[a].id].B, *(ep_t *)(ctx->dblpoint[a].P0));
            if(point52_is_infty(ctx->bucket[a].B->g1))
            {
                pair_mont64_to_mont52(ctx->bucket[a].B, *(ctx->bucket[a].wait));
            }
            else
            {
                ctx->avx3_out[ctx->avx3_cnt] = &(ctx->bucket[a].B);
                mpi_conv_64to52(u2x[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->g1->x, HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(u2y[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->g1->y, HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(v2x[ctx->avx3_cnt][0], (*(ctx->bucket[a].wait))->g2->x[0], HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(v2y[ctx->avx3_cnt][0], (*(ctx->bucket[a].wait))->g2->y[0], HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(v2x[ctx->avx3_cnt][1], (*(ctx->bucket[a].wait))->g2->x[1], HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(v2y[ctx->avx3_cnt][1], (*(ctx->bucket[a].wait))->g2->y[1], HT_NWORDS, 6); // convert to radix-52
                // mpi_conv_64to52(u2z[ctx->avx3_cnt], (*(ep_t *)(ctx->dblpoint[a].P0))->z, HT_NWORDS, 6); // convert to radix-52
                ctx->avx3_cnt++;

                if(ctx->avx3_cnt == AVX_WAY)
                {
                    simdPADDmix_G1(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
                    simdPADDmix_G2(ctx->avx3_out, v2x, v2y, ctx->avx3_cnt);
                    ctx->avx3_cnt = 0;
                }
            }
        }
    }
    if(ctx->avx3_cnt)
    {
        simdPADDmix_G1(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
        simdPADDmix_G2(ctx->avx3_out, v2x, v2y, ctx->avx3_cnt);
    }


    pair_mont52_to_mont64(ctx->bucket[BUCKETNUM - 1].B, M[0]);
    // ep_copy(M[0]->g1, ctx->bucket[BUCKETNUM - 1].B->g1);
    // ep2_copy(M[0]->g2, ctx->bucket[BUCKETNUM - 1].B->g2);
    M[0]->g1->coord = PROJC;
    M[0]->g2->coord = PROJC;
    tmp->g1->coord = PROJC;
    tmp->g2->coord = PROJC;
    ep_copy(Gu->g1, M[0]->g1);
    ep2_copy(Gu->g2, M[0]->g2);

    for (int j = 1; j < BUCKETNUM; j++) {
        pair_mont52_to_mont64(ctx->bucket[BUCKETNUM - 1 - j].B, tmp);
        ep_add_projc(M[j]->g1, tmp->g1, M[j - 1]->g1);
        ep2_add_projective(M[j]->g2, tmp->g2, M[j - 1]->g2);
        ep_add_projc(Gu->g1, Gu->g1, M[j]->g1);
        ep2_add_projective(Gu->g2, Gu->g2, M[j]->g2);
    }

    free(M);
}

void pair_window1(const bn_t *k, const pair *P, int num, pair Gu, int index, int WBITS)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

    pair_msm_ctx *ctx = (pair_msm_ctx*)malloc(sizeof(pair_msm_ctx));
    ctx->WBITS = WBITS;
    ctx->NWINS = NWINS;
    ctx->BUCKETNUM = BUCKETNUM;
    ctx->bucket =  (pair_Bucket*)malloc(BUCKETNUM * sizeof(pair_Bucket));

    pair_window(k, P, num, Gu, ctx, index);

    free(ctx->bucket);
    free(ctx);
}

void pair_finish(const pair *Gu, pair result, int NWINS, int WBITS)
{
	pair Tu[NWINS];

	ep_copy(Tu[0]->g1, Gu[NWINS - 1]->g1);
    ep2_copy(Tu[0]->g2, Gu[NWINS - 1]->g2);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1]->g1, Tu[i - 1]->g1);
            ep2_dbl_projective(Tu[i - 1]->g2, Tu[i - 1]->g2);
		}
		ep_add_projc(Tu[i]->g1, Gu[NWINS - 1 - i]->g1, Tu[i - 1]->g1);
        ep2_add_projective(Tu[i]->g2, Gu[NWINS - 1 - i]->g2, Tu[i - 1]->g2);
	}
	ep_copy(result->g1, Tu[NWINS - 1]->g1);
    ep2_copy(result->g2, Tu[NWINS - 1]->g2);
    ep2_norm_projc(result->g2, result->g2);
}

void pair_finish_old(const pair *Gu, pair result, int NWINS, int WBITS)
{
	pair Tu[NWINS];

	ep_copy(Tu[0]->g1, Gu[NWINS - 1]->g1);
    ep2_copy(Tu[0]->g2, Gu[NWINS - 1]->g2);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1]->g1, Tu[i - 1]->g1);
            ep2_dbl_projc(Tu[i - 1]->g2, Tu[i - 1]->g2);
		}
		ep_add_projc(Tu[i]->g1, Gu[NWINS - 1 - i]->g1, Tu[i - 1]->g1);
        ep2_add_projc(Tu[i]->g2, Gu[NWINS - 1 - i]->g2, Tu[i - 1]->g2);
	}
	ep_copy(result->g1, Tu[NWINS - 1]->g1);
    ep2_copy(result->g2, Tu[NWINS - 1]->g2);
}

void pair_ifma_threads(const bn_t *k, const pair *P, int num, pair result, int WBITS)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    
    pair Gu[NWINS];

    clock_t begin, end; 
    double time_used;
    begin = clock();
    #pragma omp parallel for num_threads(16) schedule(dynamic)
    for (int i = 0; i < NWINS; i=i+1)
    {
        pair_window1(k, P, num, Gu[i], i, WBITS);
    }

    pair_finish(Gu, result, NWINS, WBITS);

}


void old_pair_window(int m, pair Gu, const bn_t *k, const pair *P, int index, int WBITS)
{
	int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

	pair *B = (pair*)malloc(BUCKETNUM * sizeof(pair));
	pair *M = (pair*)malloc(BUCKETNUM * sizeof(pair));

	uint32_t b;

    for (int j = 0; j < BUCKETNUM; j++) {
        ep_set_infty(B[j]->g1);
        ep2_set_infty(B[j]->g2);
    }
    
    for (int j = 0; j < m; j++) {
        b = get_wval(k[j], index, WBITS);
        if (b != 0) {
            ep_add_projc(B[(uint32_t)b - 1]->g1, B[(uint32_t)b - 1]->g1, P[j]->g1);
            ep2_add_projc(B[(uint32_t)b - 1]->g2, B[(uint32_t)b - 1]->g2, P[j]->g2);
        }
    }
    ep_copy(M[0]->g1, B[BUCKETNUM - 1]->g1);
    ep2_copy(M[0]->g2, B[BUCKETNUM - 1]->g2);
    ep_copy(Gu->g1, M[0]->g1);
    ep2_copy(Gu->g2, M[0]->g2);

    for (int j = 1; j < BUCKETNUM; j++) {
        ep_add_projc(M[j]->g1, B[BUCKETNUM - 1 - j]->g1, M[j - 1]->g1);
        ep2_add_projc(M[j]->g2, B[BUCKETNUM - 1 - j]->g2, M[j - 1]->g2);
        ep_add_projc(Gu->g1, Gu->g1, M[j]->g1);
        ep2_add_projc(Gu->g2, Gu->g2, M[j]->g2);
    }
    free(B);
    free(M);

}

void old_pair_threads(const bn_t *k, const pair *P, int num, pair result, int WBITS)
{
    // printf("WBITS = %d\n", WBITS);
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    
    pair Gu[NWINS];
    #pragma omp parallel for num_threads(16) schedule(dynamic)
    for (int i = 0; i < NWINS; i=i+1)
    {
        old_pair_window(num, Gu[i], k, P, i, WBITS);
    }
    pair_finish_old(Gu, result, NWINS, WBITS);

}
