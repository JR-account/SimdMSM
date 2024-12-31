#include "pippenger.h"
#include "types.h"
#include <time.h>
#include <math.h>

int main()
{
    // init Relic 
    core_init();

    // init curve
    fp_prime_init();
    ep_curve_init();

    fp_param_set(B12_381);
    ep_param_set(B12_P381);

    // result G
    ep_t G;
    bn_t N;
    ep_null(G);
    ep_new(G);
    ep_set_infty(G);
    ep_curve_get_ord(N);

    int num = 1<<20;
    // int round = num / ROUNDNUM + 1;

    /*------------------init P and k-------------------------*/
    ep_t *point = (ep_t*)malloc(num * sizeof(ep_t));
    for (int i = 0; i < num; i++)
    {
        ep_rand(point[i]);
        ep_norm(point[i], point[i]);
    }

    bn_t *k1 = (bn_t*)malloc(num * sizeof(bn_t));
    for (int i = 0; i < num; i++) 
    {
        bn_null(k1[i]);
        bn_rand_mod(k1[i], N);
    }

    pip_ctx ctx;
    int pracs[] = {11, 12, 12, 13, 13, 14, 14, 14, 14, 15};
    int log2_length = log(num)/log(2);
    if(log2_length < 8)
    {
        ctx.WBITS = 5;
    } 
    else if(log2_length < 11)
    {
        ctx.WBITS = log2_length - 3;
    }
    else if(log2_length < 15)
    {
        ctx.WBITS = log2_length - 4;
    }
    else if(log2_length < 25)
    {
        ctx.WBITS = pracs[log2_length - 15];
    }
    else
    {
        ctx.WBITS = log2_length - 9;
    }
    printf("WBITS = %d\n", ctx.WBITS);
    ctx.NWINS = ((NBITS + ctx.WBITS - 1) / ctx.WBITS);
    ctx.BUCKETNUM = (1 << ctx.WBITS) - 1;
    ctx.bucket = (Bucket*)malloc(ctx.BUCKETNUM * sizeof(Bucket));


    clock_t begin, end; 
    double time_used;
    printf("Pippenger_old: ");
    begin = clock();
    for (int j = 0; j < 1; j++)
    {
        Pippenger_old(num, G, k1,  point, &ctx);
    }
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC/100;
    printf("time_used = %lf\n", time_used);
    ep_norm(G, G);
    ep_print(G);
    
    ep_t *M = (ep_t*)malloc(ctx.BUCKETNUM * sizeof(ep_t));

    init(&ctx);
    printf("Pippenger_ifma: ");
    begin = clock();
    for (int j = 0; j < 1; j++)
    {
        pippenger_fast(k1, point, num, G,  M, &ctx);
    }
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC/100;
    printf("time_used = %lf\n", time_used);


    if(ep_on_curve(G) != 1)
    {
        printf("not on curve!\n");
    }
    ep_norm(G, G);
    ep_print(G);

    free(M);
    free(ctx.bucket);
}