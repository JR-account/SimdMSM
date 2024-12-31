#include "pair_threads.h"
#include "types.h"
#include <time.h>
#include <math.h>
#include <sys/time.h>

#ifndef WNUM
// #define WNUM 15
#endif
#ifndef WMBITS
// #define WMBITS 3
#endif

#ifndef NBENCHS 
// #define NBENCHS 1
#endif

uint64_t clockMS()
{
  struct timeval te;
  gettimeofday(&te, NULL);                                                                // get current time
  uint64_t milliseconds = ((uint64_t)te.tv_sec) * 1000L + ((uint64_t)te.tv_usec) / 1000L; // calculate milliseconds
  return milliseconds;
}

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
    // int round = num / ROUNDNUM + 1;

    /*------------------init P and k-------------------------*/
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



    // FILE *pointFile = fopen("pair.txt", "r");
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

    uint64_t timer;
    double time_used;
    printf("threads_old");
    timer = -clockMS();
    for (int j = 0; j < 1; j++)
    {
        old_pair_threads(k1, point, num, G1, WMBITS);
    }
    timer += clockMS();
    time_used = (double)timer / 1000;
    printf("=%lf\n", time_used);
    ep_norm(G1->g1, G1->g1);
    // ep_print(G1->g1);
    ep2_norm(G1->g2, G1->g2);
    // ep2_print(G1->g2);
    
    printf("threads_ifma");
    timer = -clockMS();
    for (int j = 0; j < NBENCHS; j++)
    {
        pair_ifma_threads(k1, point, num, G2, WMBITS);
    }
    timer += clockMS();
    time_used = (double)timer / 1000 / NBENCHS;
    printf("=%lf\n", time_used);


    ep_norm(G2->g1, G2->g1);
    // ep_print(G2->g1);
    ep2_norm(G2->g2, G2->g2);
    // ep2_print(G2->g2);
    // if(ep_on_curve(G2->g1))
    // {
    //     printf("g1 on curve\n");
    // }
    // else
    // {
    //     printf("g1 not on curve\n");
    // }

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

}