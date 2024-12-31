#include "pippenger.h"


/* "get_wval" function is used to get sub-scalar from the curve */
uint32_t get_wval(const bn_t d, uint32_t offset, int WBITS)
{
	bn_t rs;
    bn_rsh(rs, d, offset * WBITS);
	return rs->dp[0] & ((1 << WBITS) - 1);
}

void Pippenger_old(int m, ep_t G, const bn_t *k, const ep_t *P,  pip_ctx *ctx)
{
	int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

	// ep_t B[BUCKETNUM], M[BUCKETNUM];
	ep_t *B = (ep_t*)malloc(BUCKETNUM * sizeof(ep_t));
	ep_t *M = (ep_t*)malloc(BUCKETNUM * sizeof(ep_t));

	ep_t Gu[NWINS], Tu[NWINS];
	uint32_t b;

	for (int i = 0; i < NWINS; i++) {
		for (int j = 0; j < BUCKETNUM; j++) {
			ep_set_infty(B[j]);
		}
        
		for (int j = 0; j < m; j++) {
            b = get_wval(k[j], i, WBITS);
			if (b != 0) {
				ep_add_projc(B[(uint32_t)b - 1], B[(uint32_t)b - 1], P[j]);
			}
		}
		ep_copy(M[0], B[BUCKETNUM - 1]);
		ep_copy(Gu[i], M[0]);

		for (int j = 1; j < BUCKETNUM; j++) {
			ep_add_projc(M[j], B[BUCKETNUM - 1 - j], M[j - 1]);
			ep_add_projc(Gu[i], Gu[i], M[j]);
		}
        // printf("第%d列 \n", i);
		// ep_norm(Gu[i], Gu[i]);
        // ep_print(Gu[i]);

	}
	ep_copy(Tu[0], Gu[NWINS - 1]);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1], Tu[i - 1]);
		}
		ep_add_projc(Tu[i], Gu[NWINS - 1 - i], Tu[i - 1]);
	}
	ep_copy(G, Tu[NWINS - 1]);

}


void MSM_old(int m, ep_t G, const bn_t *k, const ep_t *P)
{
    ep_t tmp;
    ep_null(tmp);
    ep_new(tmp);

    ep_set_infty(G); // 初始化结果为无穷远点

    for (int i = 0; i < m; i++)
    {
        // tmp = k[i] * P[i]
        ep_mul_basic(tmp, P[i], k[i]);

        // G = G + tmp
        ep_add_basic(G, G, tmp);

        
        // ep_free(tmp);
    }
}

void point52_set_infty(point52_t p) {
	for(int i = 0; i < HT_NWORDS; i++)
    {
        p->x[i] = 0;
        p->y[i] = 0;
        p->z[i] = 0;
    }
}

void init(pip_ctx *ctx)
{
	int i;
	/* init Bucket*/
	for(i = 0; i < ctx->BUCKETNUM; i++)
	{
		ctx->bucket[i].id = i;
		point52_set_infty(ctx->bucket[i].B);
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
		point52_set_infty(ctx->T.pbuf[i].point);
	}

	/* init avx3*/
	ctx->avx3_cnt = 0;
	for(i = 0; i < AVX_WAY; i++)
	{
		ctx->avx3_out[i] = NULL;
	}
	
}

/* avx2 is simdPADDprj*/
void avx2(pip_ctx *ctx, int eff_times)
{
	int i;
	buf2_st *item;
	ep_t tmp;

	/* once simdPADDprj*/
	for(i = 0; i < eff_times; i++)
	{
		item = &(ctx->buf2[i]);
		if(item->flag == 0)
		{
			ep_add_projc(tmp, *item->op1, *item->op2);
			queue_in(&(ctx->T), tmp, item->id);
		}
		else
		{
			ep_add_projc(*item->op2, *item->op1, *item->op2);
		}
	}
	ctx->cnt2 = 0;
}

buf2_st* check_eq_id(pip_ctx *ctx, int id) {
    int i = 0;
    while (i < ctx->cnt2) {
        if ((ctx->buf2[i].id == id) && (ctx->buf2[i].flag == 1)) {
            return &(ctx->buf2[i]);
        }
        i++;
    }
    return NULL;
}

/* process of T*/
void state_trans_T(pip_ctx *ctx)
{
	buf2_st *same_id;
	int bucket_id, i;

	while (!queue_is_empty(&(ctx->T))) 
	{
		i = ctx->T.front;
        bucket_id = ctx->T.pbuf[i].id;
		same_id = check_eq_id(ctx, bucket_id);
		if(same_id == NULL)
		{	
			/* write into bucket*/
			ctx->buf2[ctx->cnt2].id = bucket_id;
			ctx->buf2[ctx->cnt2].op1 = &(ctx->T.pbuf[i].point);
			ctx->buf2[ctx->cnt2].op2 = &(ctx->bucket[bucket_id].B);
			ctx->buf2[ctx->cnt2].flag = 1;
			// ctx->buf2[ctx->cnt2].out = &(ctx->bucket[ctx->T.pbuf[i].id].B);
			ctx->cnt2++;
		}
		else
		{
			/* write into queue*/
	    	same_id->op2 = &(ctx->T.pbuf[i].point);
			same_id->flag = 0;
			// ctx->buf2[ctx->cnt2].out = &(ctx->T.pbuf[].point);
		}

		if(ctx->cnt2 == AVX_WAY)
		{
			avx2(ctx, AVX_WAY);
		}

		queue_out(&(ctx->T), 1);
    }
}

/* avx1 is simdPADDaff*/
void avx1(pip_ctx *ctx, int eff_times)
{
	int i;
	buf1_st *dblp;
	ep_t tmp;
	/* once simdPADDaff*/
	for(i = 0; i < eff_times; i++)
	{
		dblp = &ctx->buf1[i];
		ep_add_projc(tmp, *(dblp->P0), *(dblp->P1));
		queue_in(&(ctx->T), tmp, dblp->id);
	}
	ctx->cnt1 = 0;

	/* deal with the avx1_out, change state in T*/
	state_trans_T(ctx);
	
}

void put_into_bucket(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx)
{
	/* judge if the bucket wait is empty*/
	if(ctx->bucket[b_id].wait == NULL)
	{
		// printf("%d ", p_id);
		ctx->bucket[b_id].wait = &(P[p_id]);
	}
	else
	{
		ctx->buf1[ctx->cnt1].id = b_id;
		ctx->buf1[ctx->cnt1].P0 = ctx->bucket[b_id].wait;
		ctx->buf1[ctx->cnt1].P1 = &(P[p_id]);
		ctx->cnt1++;
		ctx->bucket[b_id].wait = NULL;

		/* judge if buf1 is full*/
		if(ctx->cnt1 == AVX_WAY)
		{
			/* avx1_in is full, then do avx1*/
			avx1(ctx, AVX_WAY);
		}
	}
}


void pippenger_fast(const bn_t *k, const ep_t *P, int num, ep_t result, ep_t *M, pip_ctx *ctx)
{
	int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

	ep_t Gu[NWINS], Tu[NWINS];
	uint32_t b;
	int i, a;
	for (i = 0; i < NWINS; i++) {
		init(ctx);
		for (int j = 0; j < num; j++) {
            b = get_wval(k[j], i, WBITS);
			if (b != 0) {
				// ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
				put_into_bucket(P, (uint32_t)b - 1, j, ctx);
			}
		}

		/* deal with the remains in buf*/
		/* remain buf1*/
		if(ctx->cnt1 != 0)
		{
			avx1(ctx, ctx->cnt1);
		}

		/* remain buf2*/
		while((!queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
		{
			avx2(ctx, ctx->cnt2);
			state_trans_T(ctx);
		}


		/* remain buf0*/
		for(a = 0; a < BUCKETNUM; a++)
		{
			if(ctx->bucket[a].wait != NULL)
			{
				ep_add_projc(ctx->bucket[a].B, ctx->bucket[a].B, *(ctx->bucket[a].wait));
			}
		}

		ep_copy(M[0], ctx->bucket[BUCKETNUM - 1].B);
		ep_copy(Gu[i], M[0]);

		for (int j = 1; j < BUCKETNUM; j++) {
			ep_add_projc(M[j], ctx->bucket[BUCKETNUM - 1 - j].B, M[j - 1]);
			ep_add_projc(Gu[i], Gu[i], M[j]);
		}

	}

	ep_copy(Tu[0], Gu[NWINS - 1]);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1], Tu[i - 1]);
		}
		ep_add_projc(Tu[i], Gu[NWINS - 1 - i], Tu[i - 1]);
	}
	ep_copy(result, Tu[NWINS - 1]);
}


