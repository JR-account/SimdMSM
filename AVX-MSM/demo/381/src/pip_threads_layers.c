#include "pip_threads.h"
#include <math.h>
#include <omp.h>
#include <time.h>
#include <pthread.h>

void init2(pip_ctx *ctx)
{
	int i;
	// /* init Bucket*/
	// for(i = 0; i < ctx->BUCKETNUM; i++)
	// {
	// 	ctx->bucket[i].id = i;
	// 	point52_set_infty(ctx->bucket[i].B);
	// 	ctx->bucket[i].wait = NULL;
	// }

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

void avx2_withIFMA_threads(pip_ctx *ctx, int eff_times, pthread_mutex_t* mutex)
{
    int i;
    buf2_st *item;
    point52_t avx2_out[AVX_WAY];
    
    pthread_mutex_lock(mutex);
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
    pthread_mutex_unlock(mutex);

    ctx->cnt2 = 0;
    
}

/* process of T*/
void state_trans_T_withIFMA_threads(pip_ctx *ctx, pthread_mutex_t* mutex)
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
            // pthread_mutex_lock(mutex);
			avx2_withIFMA_threads(ctx, AVX_WAY, mutex);
            // pthread_mutex_unlock(mutex);
		}

		queue_out(&(ctx->T), 1);
    }
}


void avx1_withIFMA_threads(pip_ctx *ctx, int eff_times, pthread_mutex_t* mutex)
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
        pthread_mutex_lock(mutex);
        dblp = &ctx->buf1[i];
        if(point52_is_infty(ctx->bucket[dblp->id].B))
        {
            point52_copy(ctx->bucket[dblp->id].B, avx1_out[i]);
        }
        else
        {
            queue_in(&(ctx->T), avx1_out[i], dblp->id);
        }
        pthread_mutex_unlock(mutex);
    }

	ctx->cnt1 = 0;

	/* deal with the avx1_out, change state in T*/
	state_trans_T_withIFMA_threads(ctx, mutex);
}

void put_into_bucket_threads(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx, pthread_mutex_t* mutex)
{
	/* judge if the bucket wait is empty*/
	if(ctx->bucket[b_id].wait == NULL)
	{
		// printf("%d ", p_id);
        pthread_mutex_lock(mutex);
		ctx->bucket[b_id].wait = &(P[p_id]);
        pthread_mutex_unlock(mutex);
	}
	else
	{
        pthread_mutex_lock(mutex);
		ctx->buf1[ctx->cnt1].id = b_id;
		ctx->buf1[ctx->cnt1].P0 = ctx->bucket[b_id].wait;
		ctx->buf1[ctx->cnt1].P1 = &(P[p_id]);
		ctx->cnt1++;
		ctx->bucket[b_id].wait = NULL;
        pthread_mutex_unlock(mutex);

		/* judge if buf1 is full*/
		if(ctx->cnt1 == AVX_WAY)
		{
			/* avx1_in is full, then do avx1*/
			avx1_withIFMA_threads(ctx, AVX_WAY, mutex);
		}
	}
}

void loop_msm(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index, int loop, pthread_mutex_t* mutex, int WBITS, Bucket* bucket)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

    uint32_t b;
    pip_ctx *ctx = (pip_ctx*)malloc(sizeof(pip_ctx));
    ctx->bucket = bucket;

    init2(ctx);

    for(int j = loop*(num/2); j < (loop+1)*(num/2); j++)
    {
        b = get_wval(k[j], index, WBITS);
        if (b != 0)
        {
            // ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
            // pthread_mutex_lock(mutex);
            put_into_bucket_threads(P, (uint32_t)b - 1, j, ctx, mutex);
            // pthread_mutex_unlock(mutex);
        }
    }

        /* deal with the remains in buf*/
    /* remain buf1*/
    if(ctx->cnt1 != 0)
    {
        // pthread_mutex_lock(mutex);
        avx1_withIFMA_threads(ctx, ctx->cnt1, mutex);
        // pthread_mutex_unlock(mutex);
    }

    /* remain buf2*/
    while((!queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
    {
        // pthread_mutex_lock(mutex);
        avx2_withIFMA_threads(ctx, ctx->cnt2, mutex);
        state_trans_T_withIFMA_threads(ctx, mutex);
        // pthread_mutex_unlock(mutex);
    }
}

void submsm(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index, int WBITS, Bucket* bucket)
{
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;
    uint32_t b;
    int i, a;
    ep_t tmp;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    ep_t *M = (ep_t*)malloc(BUCKETNUM * sizeof(ep_t));
    point52_t *avx3_out[AVX_WAY];
    int avx3_cnt = 0;
    for(i = 0; i < AVX_WAY; i++)
	{
		avx3_out[i] = NULL;
	}

    // printf("第%d列 \n", i);
    pthread_mutex_t mutex;

    pthread_mutex_init(&mutex,NULL);


    #pragma omp parallel 
    for (int j = 0; j < 2; j++)
    {

        loop_msm(k, P, num, Gu, index, j, &mutex, WBITS, bucket);
    }


    /* deal with the remains in buf0*/
    for (a = 0; a < BUCKETNUM; a++)
    {
        if (bucket[a].wait != NULL)
        {
            // ep_add_projc(ctx->bucket[ctx->dblpoint[a].id].B, ctx->bucket[ctx->dblpoint[a].id].B, *(ep_t *)(ctx->dblpoint[a].P0));
            if(point52_is_infty(bucket[a].B))
            {
                mont64_to_mont52(bucket[a].B, *(bucket[a].wait));
            }
            else
            {
                avx3_out[avx3_cnt] = &(bucket[a].B);
                mpi_conv_64to52(u2x[avx3_cnt], (*(bucket[a].wait))->x, HT_NWORDS, 6); // convert to radix-52
                mpi_conv_64to52(u2y[avx3_cnt], (*(bucket[a].wait))->y, HT_NWORDS, 6); // convert to radix-52
                // mpi_conv_64to52(u2z[ctx->avx3_cnt], (*(ep_t *)(ctx->dblpoint[a].P0))->z, HT_NWORDS, 6); // convert to radix-52
                avx3_cnt++;

                if(avx3_cnt == AVX_WAY)
                {
                    simdPADDmix(avx3_out, u2x, u2y, avx3_cnt);
                    avx3_cnt = 0;
                }
            }
        }
    }
    if(avx3_cnt)
    {
        simdPADDmix(avx3_out, u2x, u2y, avx3_cnt);
    }

    /* Summary points in buckets*/
    mont52_to_mont64(bucket[BUCKETNUM - 1].B, M[0]);
    // ep_copy(M[0], ctx->bucket[BUCKETNUM - 1].B);
    M[0]->coord = PROJC;
    tmp->coord = PROJC;
    ep_copy(Gu, M[0]);

    for (int j = 1; j < BUCKETNUM; j++)
    {
        mont52_to_mont64(bucket[BUCKETNUM - 1 - j].B, tmp);
        ep_add_projc(M[j], tmp, M[j - 1]);
        ep_add_projc(Gu, Gu, M[j]);
    }
    free(M);

}

void pip_windows(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index)
{
    int pracs[] = {11, 12, 12, 13, 13, 14, 14, 14, 14, 15};
    int log2_length = log(num)/log(2);
    int WBITS = 0;
    if(log2_length < 8)
    {
        WBITS = log2_length - 2;
    } 
    else if(log2_length < 11)
    {
        WBITS = log2_length - 3;
    }
    else if(log2_length < 15)
    {
        WBITS = log2_length - 4;
    }
    else if(log2_length < 25)
    {
        WBITS = pracs[log2_length - 15];
    }
    else
    {
        WBITS = log2_length - 9;
    }
    // printf("WBITS = %d\n", WBITS);
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    int BUCKETNUM = (1 << WBITS) - 1;

    // pip_ctx *ctx = (pip_ctx*)malloc(sizeof(pip_ctx));
    // ctx->WBITS = WBITS;
    // ctx->NWINS = NWINS;
    // ctx->BUCKETNUM = BUCKETNUM;
    Bucket* bkt = (Bucket*)malloc(BUCKETNUM * sizeof(Bucket));
	/* init Bucket*/
	for(int i = 0; i < BUCKETNUM; i++)
	{
		bkt[i].id = i;
		point52_set_infty(bkt[i].B);
		bkt[i].wait = NULL;
	}
    submsm(k, P, num, Gu, index, WBITS, bkt);

    // free(ctx->bucket);
    // free(ctx);
}


void pip_threads_layers(const bn_t *k, const ep_t *P, int num, ep_t result)
{
    int pracs[] = {11, 12, 12, 13, 13, 14, 14, 14, 14, 15};
    int log2_length = log(num)/log(2);
    int WBITS = 0;
    if(log2_length < 8)
    {
        WBITS = log2_length - 2;
    } 
    else if(log2_length < 11)
    {
        WBITS = log2_length - 3;
    }
    else if(log2_length < 15)
    {
        WBITS = log2_length - 4;
    }
    else if(log2_length < 25)
    {
        WBITS = pracs[log2_length - 15];
    }
    else
    {
        WBITS = log2_length - 9;
    }
    // printf("WBITS = %d\n", WBITS);
    int NWINS = ((NBITS + WBITS - 1) / WBITS);
    
    ep_t Gu[NWINS];

    // pip_ctx ctx[NWINS];

    // for (i = 0; i < NWINS; i++)
    // {
    //     ctx[i].WBITS = WBITS;

    //     ctx[i].NWINS = NWINS;
    //     ctx[i].BUCKETNUM = BUCKETNUM;
    //     ctx[i].bucket = (Bucket*)malloc(BUCKETNUM * sizeof(Bucket));

    // }
    // clock_t begin, end; 
    // double time_used;
    // begin = clock();
    // double total_start_time = omp_get_wtime();
    #pragma omp parallel for num_threads(18) schedule(dynamic)
    for (int i = 0; i < NWINS; i=i+1)
    {
        // int thread_id = omp_get_thread_num(); // 获取线程 ID
        // double start_time = omp_get_wtime();  // 获取开始时间

        // 执行你的函数
        pip_windows(k, P, num, Gu[i], i);
        // if(i+1 < NWINS)
        // {
        //     pip_window1(k, P, num, Gu[i+1], i+1);
        // }
        // if(i+2 < NWINS)
        // {
        //     pip_window1(k, P, num, Gu[i+2], i+2);
        // }


        // pip_window1(k, P, num, Gu[i+1], i+1);

        // double end_time = omp_get_wtime(); // 获取结束时间
        // double elapsed_time = (end_time - start_time);

        // printf("Thread ID: %d, Execution Time: %.16g seconds\n", thread_id, elapsed_time);
    }
    // // 记录并行区域结束时间
    // double total_end_time = omp_get_wtime();
    // double total_elapsed_time = total_end_time - total_start_time;

    // // 输出整个并行区域的总时间
    // printf("Total Execution Time for Parallel Region: %.16g seconds\n", total_elapsed_time);

    // // for (int i = 9; i < NWINS; i++)
    // {
    //     // pip_ctx ctx;
    //     // ctx.WBITS = WBITS;
    //     // ctx.NWINS = NWINS;
    //     // ctx.BUCKETNUM = BUCKETNUM;
    //     // ctx.bucket = (Bucket*)malloc(BUCKETNUM * sizeof(Bucket));
    //     // pip_window(k, P, num, Gu[i], &ctx, i);
    //     pip_window1(k, P, num, Gu[i], i);
    // }

    // #pragma omp parallel for num_threads(3)
    // for (i = 12; i < NWINS; i++)
    // {
    //     pip_window(k, P, num, Gu[i], &ctx[i], i);
    // }

    // finish(Gu, result, NWINS, WBITS);

    // clock_t begin, end; 
    // double time_used;
    // begin = clock();

    // // 记录 finish 函数开始时间
    // double finish_start_time = omp_get_wtime();

    // 调用 finish 函数
    finish(Gu, result, NWINS, WBITS);

    // // 记录 finish 函数结束时间
    // double finish_end_time = omp_get_wtime();
    // double finish_elapsed_time = finish_end_time - finish_start_time;

    // // 输出 finish 函数的执行时间
    // printf("Execution Time for finish: %.16g seconds\n", finish_elapsed_time);

    // end = clock();
    // time_used = (double)(end - begin) / CLOCKS_PER_SEC;
    // printf("time_used = %lf\n", time_used);

}
