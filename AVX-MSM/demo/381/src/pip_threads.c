#include "pip_threads.h"
#include <math.h>
#include <omp.h>
#include <time.h>

void pip_window(const bn_t *k, const ep_t *P, int num, ep_t Gu, pip_ctx *ctx, int index)
{
    int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;
    uint32_t b;
    int i, a;
    ep_t tmp;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    ep_t *M = (ep_t*)malloc(ctx->BUCKETNUM * sizeof(ep_t));

    // printf("第%d列 \n", i);
    init(ctx);
    for (int j = 0; j < num; j++)
    {
        b = get_wval(k[j], index, ctx->WBITS);
        if (b != 0)
        {
            // ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
            put_into_bucket_withIFMA(P, (uint32_t)b - 1, j, ctx);
        }
    }

    /* deal with the remains in buf*/
    /* remain buf1*/
    if(ctx->cnt1 != 0)
    {
        avx1_withIFMA(ctx, ctx->cnt1);
    }

    /* remain buf2*/
    while((!queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
    {
        avx2_withIFMA(ctx, ctx->cnt2);
        state_trans_T_withIFMA(ctx);
    }

    /* deal with the remains in buf0*/
    for (a = 0; a < BUCKETNUM; a++)
    {
        if (ctx->bucket[a].wait != NULL)
        {
            // ep_add_projc(ctx->bucket[ctx->dblpoint[a].id].B, ctx->bucket[ctx->dblpoint[a].id].B, *(ep_t *)(ctx->dblpoint[a].P0));
            if(point52_is_infty(ctx->bucket[a].B))
            {
                mont64_to_mont52(ctx->bucket[a].B, *(ctx->bucket[a].wait));
            }
            else
            {
                ctx->avx3_out[ctx->avx3_cnt] = &(ctx->bucket[a].B);
                mpi_conv_64to52(u2x[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->x, HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(u2y[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->y, HT_NWORDS, 6); // convert to radix-52
                // mpi_conv_64to52(u2z[ctx->avx3_cnt], (*(ep_t *)(ctx->dblpoint[a].P0))->z, HT_NWORDS, 6); // convert to radix-52
                ctx->avx3_cnt++;

                if(ctx->avx3_cnt == AVX_WAY)
                {
                    simdPADDmix(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
                    ctx->avx3_cnt = 0;
                }
            }
        }
    }
    if(ctx->avx3_cnt)
    {
        simdPADDmix(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
    }

    /* Summary points in buckets*/
    mont52_to_mont64(ctx->bucket[BUCKETNUM - 1].B, M[0]);
    // ep_copy(M[0], ctx->bucket[BUCKETNUM - 1].B);
    M[0]->coord = PROJC;
    tmp->coord = PROJC;
    ep_copy(Gu, M[0]);

    for (int j = 1; j < BUCKETNUM; j++)
    {
        mont52_to_mont64(ctx->bucket[BUCKETNUM - 1 - j].B, tmp);
        ep_add_projc(M[j], tmp, M[j - 1]);
        ep_add_projc(Gu, Gu, M[j]);
    }
    free(M);

}

void pip_window1(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index, int WBITS)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

    pip_ctx *ctx = (pip_ctx*)malloc(sizeof(pip_ctx));
    ctx->WBITS = WBITS;
    ctx->NWINS = NWINS;
    ctx->BUCKETNUM = BUCKETNUM;
    ctx->bucket = (Bucket*)malloc(BUCKETNUM * sizeof(Bucket));

    pip_window(k, P, num, Gu, ctx, index);

    free(ctx->bucket);
    free(ctx);
}

void finish(const ep_t *Gu, ep_t result, int NWINS, int WBITS)
{
	ep_t Tu[NWINS];

	ep_copy(Tu[0], Gu[NWINS - 1]);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1], Tu[i - 1]);
		}
		ep_add_projc(Tu[i], Gu[NWINS - 1 - i], Tu[i - 1]);
	}
	ep_copy(result, Tu[NWINS - 1]);
}

void pip_threads(const bn_t *k, const ep_t *P, int num, ep_t result, int WBITS)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    
    ep_t Gu[NWINS];

    #pragma omp parallel for num_threads(18) schedule(dynamic)
    for (int i = 0; i < NWINS; i=i+1)
    {
        pip_window1(k, P, num, Gu[i], i, WBITS);

    }
    finish(Gu, result, NWINS, WBITS);
}


void old_window(int m, ep_t Gu, const bn_t *k, const ep_t *P, int index, int WBITS)
{
	int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

	// ep_t B[BUCKETNUM], M[BUCKETNUM];
	ep_t *B = (ep_t*)malloc(BUCKETNUM * sizeof(ep_t));
	ep_t *M = (ep_t*)malloc(BUCKETNUM * sizeof(ep_t));

	uint32_t b;


    for (int j = 0; j < BUCKETNUM; j++) {
        ep_set_infty(B[j]);
    }
    
    for (int j = 0; j < m; j++) {
        b = get_wval(k[j], index, WBITS);
        if (b != 0) {
            ep_add_projc(B[(uint32_t)b - 1], B[(uint32_t)b - 1], P[j]);
        }
    }
    ep_copy(M[0], B[BUCKETNUM - 1]);
    ep_copy(Gu, M[0]);

    for (int j = 1; j < BUCKETNUM; j++) {
        ep_add_projc(M[j], B[BUCKETNUM - 1 - j], M[j - 1]);
        ep_add_projc(Gu, Gu, M[j]);
    }
    free(B);
    free(M);


}

void old_threads(const bn_t *k, const ep_t *P, int num, ep_t result, int WBITS)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    
    ep_t Gu[NWINS];

    #pragma omp parallel for num_threads(16) schedule(dynamic)
    for (int i = 0; i < NWINS; i=i+1)
    {
        old_window(num, Gu[i], k, P, i, WBITS);
    }

    finish(Gu, result, NWINS, WBITS);

}
