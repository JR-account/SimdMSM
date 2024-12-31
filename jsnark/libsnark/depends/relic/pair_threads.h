#ifndef PAIR_THREADS_H
#define PAIR_THREADS_H

#include "pair_ifma.h"
#include <math.h>

void pair_window(const bn_t *k, const pair *P, int num, pair Gu, pair_msm_ctx *ctx, int index);
void pair_window1(const bn_t *k, const pair *P, int num, pair Gu, int index, int WBITS);
void pair_finish(const pair *Gu, pair result, int NWINS, int WBITS);
void pair_ifma_threads(const bn_t *k, const pair *P, int num, pair result, int WBITS);
void old_pair_window(int m, pair Gu, const bn_t *k, const pair *P, int index, int WBITS);
void pair_finish_old(const pair *Gu, pair result, int NWINS, int WBITS);
void old_pair_threads(const bn_t *k, const pair *P, int num, pair result, int WBITS);
#endif // PIP_THREADS_H