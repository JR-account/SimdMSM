#include "pair_msm.h"
#include "pippenger.h"


void Pippenger_pair_old(int m, pair G, const bn_t *k, const pair *P, pair_msm_ctx *ctx)
{
	int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

	// pair B[BUCKETNUM], M[BUCKETNUM];
	pair Gu[NWINS], Tu[NWINS];
	pair *B = (pair*)malloc(BUCKETNUM * sizeof(pair));
	pair *M = (pair*)malloc(BUCKETNUM * sizeof(pair));

	uint32_t b;

	for (int i = 0; i < NWINS; i++) {
		for (int j = 0; j < BUCKETNUM; j++) {
			ep_set_infty(B[j]->g1);
            ep2_set_infty(B[j]->g2);
		}
        
		for (int j = 0; j < m; j++) {
            b = get_wval(k[j], i, WBITS);
			if (b != 0) {
				ep_add_projc(B[(uint32_t)b - 1]->g1, B[(uint32_t)b - 1]->g1, P[j]->g1);
                ep2_add_projc(B[(uint32_t)b - 1]->g2, B[(uint32_t)b - 1]->g2, P[j]->g2);
			}
		}
		ep_copy(M[0]->g1, B[BUCKETNUM - 1]->g1);
        ep2_copy(M[0]->g2, B[BUCKETNUM - 1]->g2);
        ep_copy(Gu[i]->g1, M[0]->g1);
        ep2_copy(Gu[i]->g2, M[0]->g2);
		// ep_copy(Gu[i], M[0]);

		for (int j = 1; j < BUCKETNUM; j++) {
			ep_add_projc(M[j]->g1, B[BUCKETNUM - 1 - j]->g1, M[j - 1]->g1);
            ep2_add_projc(M[j]->g2, B[BUCKETNUM - 1 - j]->g2, M[j - 1]->g2);
            ep_add_projc(Gu[i]->g1, Gu[i]->g1, M[j]->g1);
            ep2_add_projc(Gu[i]->g2, Gu[i]->g2, M[j]->g2);
			// ep_add_projc(Gu[i], Gu[i], M[j]);
		}

		// printf("Gu[%d]:\n", i);
		// ep2_norm(Gu[i]->g2,Gu[i]->g2);
		// ep2_print(Gu[i]->g2);
	}
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
	ep_copy(G->g1, Tu[NWINS - 1]->g1);
    ep2_copy(G->g2, Tu[NWINS - 1]->g2);

}


void MSM_pair_old(int m, pair G, const bn_t *k, const pair *P)
{
    pair tmp;

    ep_set_infty(G->g1); 
    ep2_set_infty(G->g2); 

    for (int i = 0; i < m; i++)
    {
        // tmp = k[i] * P[i]
        ep_mul_basic(tmp->g1, P[i]->g1, k[i]);
        ep2_mul_basic(tmp->g2, P[i]->g2, k[i]);

        // G = G + tmp
        ep_add_basic(G->g1, G->g1, tmp->g1);
        ep2_add_basic(G->g2, G->g2, tmp->g2);

        // ep_free(tmp);
    }
}

void pair52_set_infty(pair52 p) {
	for(int i = 0; i < HT_NWORDS; i++)
    {
        p->g1->x[i] = 0;
        p->g1->y[i] = 0;
        p->g1->z[i] = 0;
		p->g2->x[0][i] = 0;
		p->g2->x[1][i] = 0;
		p->g2->y[0][i] = 0;
		p->g2->y[1][i] = 0;
		p->g2->z[0][i] = 0;
		p->g2->z[1][i] = 0;
    }
}


void init_pair_ctx(pair_msm_ctx *ctx)
{
    int BUCKETNUM = ctx->BUCKETNUM;
	int i;
	/* init Bucket*/
	for(i = 0; i < BUCKETNUM; i++)
	{
		ctx->bucket[i].id = i;
		pair52_set_infty(ctx->bucket[i].B);
		ctx->bucket[i].wait = NULL;
	}

	/* init buf1*/
	ctx->cnt1 = 0;
	for(i = 0; i < AVX_WAY; i++)
	{
		ctx->buf1[i].id = 0;
		ctx->buf1[i].P0 = NULL;
		ctx->buf1[i].P1 = NULL;
	}

	/* init buf2*/
	ctx->cnt2 = 0;
	for(i = 0; i < AVX_WAY; i++)
	{
		ctx->buf2[i].id = 0;
		ctx->buf2[i].op1 = NULL;
		ctx->buf2[i].op2 = NULL;
		ctx->buf2[i].flag = 1;
	}

	/* init T*/
	ctx->T.front = 0;
	ctx->T.rear = 0;
	for(i = 0; i < QUEUE_SIZE; i++)
	{
		ctx->T.pbuf[i].id = 0;
		pair52_set_infty(ctx->T.pbuf[i].point);
	}

	/* init avx3*/
	ctx->avx3_cnt = 0;
	for(i = 0; i < AVX_WAY; i++)
	{
		ctx->avx3_out[i] = NULL;
	}

}
