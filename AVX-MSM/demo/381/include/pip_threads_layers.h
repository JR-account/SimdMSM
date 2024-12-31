#ifndef PIP_THREADS_H
#define PIP_THREADS_H

#include "pip_ifma.h"
#include <math.h>

void init2(pip_ctx *ctx);
void loop_msm(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index, int loop, pthread_mutex_t mutex, int WBITS, Bucket* bucket);
void submsm(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index, int WBITS, Bucket* bucket);
void pip_windows(const bn_t *k, const ep_t *P, int num, ep_t Gu, int index);
void pip_threads_layers(const bn_t *k, const ep_t *P, int num, ep_t result);
#endif // PIP_THREADS_H