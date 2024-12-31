#ifndef PIP_THREADS_H
#define PIP_THREADS_H

#include "pip_ifma.h"
#include <math.h>

void pip_window(const bn_t *k, const ep_t *P, int num, ep_t Gu, pip_ctx *ctx, int index);
void finish(const ep_t *Gu, ep_t result, int NWINS, int WBITS);
void pip_threads(const bn_t *k, const ep_t *P, int num, ep_t result, int WBITS);
void old_window(int m, ep_t Gu, const bn_t *k, const ep_t *P, int index, int WBITS);
void old_threads(const bn_t *k, const ep_t *P, int num, ep_t result, int WBITS);
#endif // PIP_THREADS_H