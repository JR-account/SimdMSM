#include "pair_ifma.h"
#include <time.h>
#include <math.h>

#ifndef WNUM
#define WNUM 10
#endif
#ifndef WMBITS
#define WMBITS 3
#endif

#ifndef NBENCHS 
#define NBENCHS 1
#endif


int main()
{
    // 初始化 Relic 库
    core_init();

    // 初始化椭圆曲线参数
    fp_prime_init();
    ep_curve_init();
    ep2_curve_init();
    fp_param_set(B12_381);
    ep_param_set(B12_P381);
    
    int type = RLC_EP_MTYPE;
    ep2_curve_set_twist(type);
    bn_t N;
    ep_curve_get_ord(N);
    // ep2_curve_get_ord(N);
    // bn_print(N);
    
    //结果G
    pair G1, G2;
    ep_set_infty(G1->g1);
    ep2_set_infty(G1->g2);
    ep_set_infty(G2->g1);
    ep2_set_infty(G2->g2);
    

    int num = (1<<WNUM);

    /*------------------init P and k-------------------------*/
    pair *point = (pair*)malloc(num * sizeof(pair));
    for (int i = 0; i < num; i++)
    {
        ep_rand(point[i]->g1);
        ep2_rand(point[i]->g2);
    }

    bn_t *k1 = (bn_t*)malloc(num * sizeof(bn_t));

    for (int i = 0; i < num; i++) 
    {
        bn_null(k1[i]);
        bn_rand_mod(k1[i], N);
    }
    
    // pair *point = (pair*)malloc(num * sizeof(pair));
    // bn_t *k1 = (bn_t*)malloc(num * sizeof(bn_t));
    //     FILE *pointFile = fopen("pair.txt", "r");
    // FILE *scalarFile = fopen("scalar.txt", "r");
    // if (pointFile == NULL) {
    //     printf("file pair error \n");
    //     return;
    // }
    // if (scalarFile == NULL) {
    //     printf("file scalar error \n");
    //     return;
    // }

    // for(int i = 0; i < (num/16); i++)
    // {
    //     fread(point[i*16], sizeof(pair), 16, pointFile);
    //     fread(k1[i*16], sizeof(bn_t), 16, scalarFile);
    // }
    // fclose(pointFile);
    // fclose(scalarFile);

    pair_msm_ctx ctx;
    ctx.WBITS = WMBITS;
     
    ctx.NWINS = ((NBITS + ctx.WBITS - 1) / ctx.WBITS);
    ctx.BUCKETNUM = (1 << ctx.WBITS) - 1;
    ctx.bucket = (pair_Bucket*)malloc(ctx.BUCKETNUM * sizeof(pair_Bucket));

    clock_t begin, end; 
    double time_used;
    printf("Pair_old");
    begin = clock();
    for (int j = 0; j < 1; j++)
    {
        Pippenger_pair_old(num, G1, k1,  point, &ctx);
    }
    // Pippenger_pair_old(num, G, k1,  point, &ctx);
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC /1;
    printf("=%f\n", time_used);
    ep_norm(G1->g1, G1->g1);
    // ep_print(G1->g1);
    ep2_norm(G1->g2, G1->g2);
    // ep2_print(G1->g2);
    

    pair *M = (pair*)malloc(ctx.BUCKETNUM * sizeof(pair));

    init_pair_ctx(&ctx);
    printf("Pair_ifma");
    begin = clock();
    for (int j = 0; j < NBENCHS; j++)
    {
        pippenger_pair_ifma(k1, point, num, G2,  M, &ctx);
    }
    end = clock();
    time_used = (double)(end - begin) / CLOCKS_PER_SEC/ NBENCHS;
    printf("=%f\n", time_used);

    ep_norm(G2->g1, G2->g1);
    // ep_print(G2->g1);
    // if(ep_on_curve(G2->g1))
    // {
    //     printf("g1 on curve\n");
    // }
    // else
    // {
    //     printf("g1 not on curve\n");
    // }
    ep2_norm(G2->g2, G2->g2);
    // ep2_print(G2->g2);
    // if(ep2_on_curve(G2->g2))
    // {
    //     printf("g2 on curve\n");
    // }
    // else
    // {
    //     printf("g2 not on curve\n");
    // }
    if((ep_cmp(G1->g1, G2->g1) == RLC_EQ) && (ep2_cmp(G1->g2, G2->g2) == RLC_EQ))
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