#include "pip_ifma.h"

void simdPADDaff(uint64_t u1x[][HT_NWORDS], uint64_t u1y[][HT_NWORDS], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], point52_t p3[8])
{
    htpoint a, b, c;
    int i;

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(u1x[7][i], u1x[6][i], u1x[5][i], u1x[4][i], u1x[3][i], u1x[2][i], u1x[1][i], u1x[0][i]);
        a.y[i] = set_vector(u1y[7][i], u1y[6][i], u1y[5][i], u1y[4][i], u1y[3][i], u1y[2][i], u1y[1][i], u1y[0][i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
    }

    gfp_mont_relic2avx_8x1w(a.x, a.x);
    gfp_mont_relic2avx_8x1w(a.y, a.y);
    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.z[i] = VSET1(ht_montR[i]);
        b.z[i] = VSET1(ht_montR[i]);
    }
    PADD_projc_8x1w_z1eqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->x, c.x, 0);
    get_channel_8x1w(p3[1]->x, c.x, 1);
    get_channel_8x1w(p3[2]->x, c.x, 2);
    get_channel_8x1w(p3[3]->x, c.x, 3);
    get_channel_8x1w(p3[4]->x, c.x, 4);
    get_channel_8x1w(p3[5]->x, c.x, 5);
    get_channel_8x1w(p3[6]->x, c.x, 6);
    get_channel_8x1w(p3[7]->x, c.x, 7);

    get_channel_8x1w(p3[0]->y, c.y, 0);
    get_channel_8x1w(p3[1]->y, c.y, 1);
    get_channel_8x1w(p3[2]->y, c.y, 2);
    get_channel_8x1w(p3[3]->y, c.y, 3);
    get_channel_8x1w(p3[4]->y, c.y, 4);
    get_channel_8x1w(p3[5]->y, c.y, 5);
    get_channel_8x1w(p3[6]->y, c.y, 6);
    get_channel_8x1w(p3[7]->y, c.y, 7);

    get_channel_8x1w(p3[0]->z, c.z, 0);
    get_channel_8x1w(p3[1]->z, c.z, 1);
    get_channel_8x1w(p3[2]->z, c.z, 2);
    get_channel_8x1w(p3[3]->z, c.z, 3);
    get_channel_8x1w(p3[4]->z, c.z, 4);
    get_channel_8x1w(p3[5]->z, c.z, 5);
    get_channel_8x1w(p3[6]->z, c.z, 6);
    get_channel_8x1w(p3[7]->z, c.z, 7);
    
}

void simdPADDprj(buf2_st* buf, point52_t p3[8])
{
    htpoint a, b, c;
    int i;

    for(i = 1; i < HT_NWORDS; i++)
    {
        if(buf[i].op1 == NULL)
        {
            buf[i].op1 = buf[0].op1;
            buf[i].op2 = buf[0].op2;
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector((*buf[7].op1)->x[i], (*buf[6].op1)->x[i], (*buf[5].op1)->x[i], (*buf[4].op1)->x[i], (*buf[3].op1)->x[i], (*buf[2].op1)->x[i], (*buf[1].op1)->x[i], (*buf[0].op1)->x[i]);
        a.y[i] = set_vector((*buf[7].op1)->y[i], (*buf[6].op1)->y[i], (*buf[5].op1)->y[i], (*buf[4].op1)->y[i], (*buf[3].op1)->y[i], (*buf[2].op1)->y[i], (*buf[1].op1)->y[i], (*buf[0].op1)->y[i]);
        a.z[i] = set_vector((*buf[7].op1)->z[i], (*buf[6].op1)->z[i], (*buf[5].op1)->z[i], (*buf[4].op1)->z[i], (*buf[3].op1)->z[i], (*buf[2].op1)->z[i], (*buf[1].op1)->z[i], (*buf[0].op1)->z[i]);
        b.x[i] = set_vector((*buf[7].op2)->x[i], (*buf[6].op2)->x[i], (*buf[5].op2)->x[i], (*buf[4].op2)->x[i], (*buf[3].op2)->x[i], (*buf[2].op2)->x[i], (*buf[1].op2)->x[i], (*buf[0].op2)->x[i]);
        b.y[i] = set_vector((*buf[7].op2)->y[i], (*buf[6].op2)->y[i], (*buf[5].op2)->y[i], (*buf[4].op2)->y[i], (*buf[3].op2)->y[i], (*buf[2].op2)->y[i], (*buf[1].op2)->y[i], (*buf[0].op2)->y[i]);
        b.z[i] = set_vector((*buf[7].op2)->z[i], (*buf[6].op2)->z[i], (*buf[5].op2)->z[i], (*buf[4].op2)->z[i], (*buf[3].op2)->z[i], (*buf[2].op2)->z[i], (*buf[1].op2)->z[i], (*buf[0].op2)->z[i]);
    }


    PADD_projc_8x1w_z1neqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->x, c.x, 0);
    get_channel_8x1w(p3[1]->x, c.x, 1);
    get_channel_8x1w(p3[2]->x, c.x, 2);
    get_channel_8x1w(p3[3]->x, c.x, 3);
    get_channel_8x1w(p3[4]->x, c.x, 4);
    get_channel_8x1w(p3[5]->x, c.x, 5);
    get_channel_8x1w(p3[6]->x, c.x, 6);
    get_channel_8x1w(p3[7]->x, c.x, 7);

    get_channel_8x1w(p3[0]->y, c.y, 0);
    get_channel_8x1w(p3[1]->y, c.y, 1);
    get_channel_8x1w(p3[2]->y, c.y, 2);
    get_channel_8x1w(p3[3]->y, c.y, 3);
    get_channel_8x1w(p3[4]->y, c.y, 4);
    get_channel_8x1w(p3[5]->y, c.y, 5);
    get_channel_8x1w(p3[6]->y, c.y, 6);
    get_channel_8x1w(p3[7]->y, c.y, 7);

    get_channel_8x1w(p3[0]->z, c.z, 0);
    get_channel_8x1w(p3[1]->z, c.z, 1);
    get_channel_8x1w(p3[2]->z, c.z, 2);
    get_channel_8x1w(p3[3]->z, c.z, 3);
    get_channel_8x1w(p3[4]->z, c.z, 4);
    get_channel_8x1w(p3[5]->z, c.z, 5);
    get_channel_8x1w(p3[6]->z, c.z, 6);
    get_channel_8x1w(p3[7]->z, c.z, 7);
    
}

void simdPADDmix(point52 *bucket[AVX_WAY], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], int num)
{
    htpoint a, b, c;
    int i;

    for(i = num; i < HT_NWORDS; i++)
    {
        if(bucket[i] == NULL)
        {
            bucket[i] = bucket[0];
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(bucket[7]->x[i], bucket[6]->x[i], bucket[5]->x[i], bucket[4]->x[i], bucket[3]->x[i], bucket[2]->x[i], bucket[1]->x[i], bucket[0]->x[i]);
        a.y[i] = set_vector(bucket[7]->y[i], bucket[6]->y[i], bucket[5]->y[i], bucket[4]->y[i], bucket[3]->y[i], bucket[2]->y[i], bucket[1]->y[i], bucket[0]->y[i]);
        a.z[i] = set_vector(bucket[7]->z[i], bucket[6]->z[i], bucket[5]->z[i], bucket[4]->z[i], bucket[3]->z[i], bucket[2]->z[i], bucket[1]->z[i], bucket[0]->z[i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
        b.z[i] = VSET1(ht_montR[i]);
   }

    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);
    PADD_projc_8x1w_mix(&c, &a, &b);

    for (i = 0; i < num ; i++)
    {
        get_channel_8x1w(bucket[i]->x, c.x, i);
        get_channel_8x1w(bucket[i]->y, c.y, i);
        get_channel_8x1w(bucket[i]->z, c.z, i);
    }
    
}


void avx2_withIFMA(pip_ctx *ctx, int eff_times)
{
    int i;
    buf2_st *item;
    point52_t avx2_out[AVX_WAY];
    
    /* avx2*/
    simdPADDprj(ctx->buf2, avx2_out);

    for (i = 0; i < eff_times; i++)
    {
        item = &(ctx->buf2[i]);
        if(item->flag == 0)
		{
			queue_in(&(ctx->T), avx2_out[i], item->id);
		}
		else
		{
			// ep_copy(ctx->bucket[item->id].B, avx2_out[i]);
            point52_copy(*item->op2, avx2_out[i]);
		}
        
    }

    ctx->cnt2 = 0;
    
}

/* process of T*/
void state_trans_T_withIFMA(pip_ctx *ctx)
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
            ctx->cnt2++;

		}
		else
		{
			/* write into queue*/
	    	same_id->op2 = &(ctx->T.pbuf[i].point);
			same_id->flag = 0;
		}

		if(ctx->cnt2 == AVX_WAY)
		{
			avx2_withIFMA(ctx, AVX_WAY);
		}

		queue_out(&(ctx->T), 1);
    }
}

int point52_is_infty(const point52_t p) {
	return (fp_is_zero(p->z) == 1);
}



void avx1_withIFMA(pip_ctx *ctx, int eff_times)
{
	int i;
	buf1_st *dblp;
	ep_t tmp;
    ep_t *p1, *p2;
    point52_t avx1_out[AVX_WAY];
    uint64_t u1x[AVX_WAY][HT_NWORDS], u1y[AVX_WAY][HT_NWORDS], u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];

	/* once simdPADDaff*/
	for(i = 0; i < eff_times; i++)
	{
		dblp = &ctx->buf1[i];
        p1 = (ep_t *)(dblp->P0);
        p2 = (ep_t *)(dblp->P1);
        mpi_conv_64to52(u1x[i], (*p1)->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u1y[i], (*p1)->y, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2x[i], (*p2)->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2y[i], (*p2)->y, HT_NWORDS, 6); // convert to radix-52
		// ep_add_projc(tmp, *(dblp->P0), *(dblp->P1));
	}
    simdPADDaff(u1x, u1y, u2x, u2y, avx1_out);
	
    for(i = 0; i < eff_times; i++)
    {
        dblp = &ctx->buf1[i];
        if(point52_is_infty(ctx->bucket[dblp->id].B))
        {
            point52_copy(ctx->bucket[dblp->id].B, avx1_out[i]);
        }
        else
        {
            queue_in(&(ctx->T), avx1_out[i], dblp->id);
        }
    }

	ctx->cnt1 = 0;

	/* deal with the avx1_out, change state in T*/
	state_trans_T_withIFMA(ctx);
}

void put_into_bucket_withIFMA(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx)
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
			avx1_withIFMA(ctx, AVX_WAY);
		}
	}
}
void mont64_to_mont52(point52_t q, ep_t p)
{
    ep_t tmp;
    fp_t mont_a;
    mont_a[0] = 0x44f6480ea8e9b9af;
    mont_a[1] = 0xa96f7d65766c8fe4;
    mont_a[2] = 0xe82efd4228b540fe;
    mont_a[3] = 0x6723e5f0ade53b2e;
    mont_a[4] = 0x25ff6eb6fdd4230a;
    mont_a[5] = 0x14c8ee06ef23c24a;

    ep_copy(tmp, p);
    fp_mul(tmp->x, tmp->x, mont_a);
    fp_mul(tmp->y, tmp->y, mont_a);
    fp_mul(tmp->z, tmp->z, mont_a);
    mpi_conv_64to52(q->x, tmp->x, HT_NWORDS, 6);
    mpi_conv_64to52(q->y, tmp->y, HT_NWORDS, 6);
    mpi_conv_64to52(q->z, tmp->z, HT_NWORDS, 6);
}

void mont52_to_mont64(point52_t q, ep_t p)
{
    ep_t tmp;
    fp_t mont_b;

    mont_b[0] = 0x0;
    mont_b[1] = 0x0;
    mont_b[2] = 0x0;
    mont_b[3] = 0x0;
    mont_b[4] = 0x0;
    mont_b[5] = 0x100000000;
    ep_set_infty(tmp);

    mpi_conv_52to64(tmp->x, q->x, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->y, q->y, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->z, q->z, 6, HT_NWORDS);
    fp_mul(p->x, tmp->x, mont_b);
    fp_mul(p->y, tmp->y, mont_b);
    fp_mul(p->z, tmp->z, mont_b); 
    // if(!ep_is_infty(p))
    // {
    //     p->coord = 2;
    // }
    // else
    // {
    //     p->coord = BASIC;
    // }

    // if(!ep_on_curve(p))
    // {
    //     printf("on curve error\n");
    // }
}
void pippenger_ifma(const bn_t *k, const ep_t *P, int num, ep_t result, ep_t *M, pip_ctx *ctx)
{
    int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

    ep_t Gu[NWINS], Tu[NWINS];
    uint32_t b;
    int i, a;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    ep_t tmp;

    for (i = 0; i < NWINS; i++)
    {
        // printf("第%d列 \n", i);
        init(ctx);
        for (int j = 0; j < num; j++)
        {
            b = get_wval(k[j], i, ctx->WBITS);
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
        M[0]->coord = PROJC;
        ep_copy(Gu[i], M[0]);
        tmp->coord = PROJC;

        for (int j = 1; j < BUCKETNUM; j++)
        {
            mont52_to_mont64(ctx->bucket[BUCKETNUM - 1 - j].B, tmp);
            ep_add_projc(M[j], tmp, M[j - 1]);
            ep_add_projc(Gu[i], Gu[i], M[j]);
        }
        // printf("第%d列 \n", i);
        // ep_norm(Gu[i], Gu[i]);
        // ep_print(Gu[i]);
        // if(Gu[i]->coord != PROJC)
        // {
        //     printf("Gu[i]->coord = \n", Gu[i]->coord);
        // }
        // if(!ep_on_curve(Gu[i]))
        // {
        //     printf("% d on curve error\n", i);
        // }



    }

    ep_copy(Tu[0], Gu[NWINS - 1]);
    for (int i = 1; i < NWINS; i++)
    {
        for (int j = 0; j < WBITS; j++)
        {
            ep_dbl_projc(Tu[i - 1], Tu[i - 1]);
        }
        ep_add_projc(Tu[i], Gu[NWINS - 1 - i], Tu[i - 1]);
    }
    ep_copy(result, Tu[NWINS - 1]);
}