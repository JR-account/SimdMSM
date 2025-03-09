#include "pip_ifma.h"
#include "types.h"
#include <time.h>
#include <math.h>

#ifndef WNUM
#define WNUM 16
#endif
#ifndef WMBITS
#define WMBITS 10
#endif

#ifndef NBENCHS 
#define NBENCHS 1
#endif

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
    ep_t G1, G2;
    bn_t N;
    ep_null(G1);
    ep_new(G1);
    ep_set_infty(G1);
    ep_null(G2);
    ep_new(G2);
    ep_set_infty(G2);
    ep_curve_get_ord(N);

    int num = (1<<WNUM);
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
    
    // FILE *pointFile = fopen("point.txt", "r");
    // FILE *scalarFile = fopen("scalar.txt", "r");
    // if (pointFile == NULL || scalarFile == NULL) {
    //     printf("file error \n");
    //     return;
    // }
    
    // for(int i = 0; i < (num/64); i++)
    // {
    //     fread(point[i*64], sizeof(ep_t), 64, pointFile);
    //     fread(k1[i*64], sizeof(bn_t), 64, scalarFile);
    // }
    // fclose(pointFile);
    // fclose(scalarFile);

    pip_ctx ctx;
    ctx.WBITS = WMBITS;

    ctx.NWINS = ((NBITS + ctx.WBITS - 1) / ctx.WBITS);
    ctx.BUCKETNUM = (1 << ctx.WBITS) - 1;
    ctx.bucket = (Bucket*)malloc(ctx.BUCKETNUM * sizeof(Bucket));

    clock_t begin, end; 
    double time_used;
    printf("Pippenger_old");
    begin = clock();
    for (int j = 0; j < 1; j++)
    {
        Pippenger_old(num, G1, k1,  point, &ctx);
    }
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("=%lf\n", time_used);
    ep_norm(G1, G1);
    // ep_print(G1);
    
    ep_t *M = (ep_t*)malloc(ctx.BUCKETNUM * sizeof(ep_t));

    init(&ctx);
    printf("Pippenger_ifma");
    begin = clock();
    for (int j = 0; j < NBENCHS; j++)
    {
        pippenger_ifma(k1, point, num, G2,  M, &ctx);
    }
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC/NBENCHS;
    printf("=%lf\n", time_used);

    // if(ep_on_curve(G2) != 1)
    // {
    //     printf("not on curve!\n");
    // }
    ep_norm(G2, G2);
    // ep_print(G2);

    if(ep_cmp(G1, G2) == RLC_EQ)
    {
        printf("YES");
    }
    else
    {
        printf("NO");
    }
    free(M);
    free(ctx.bucket);
}