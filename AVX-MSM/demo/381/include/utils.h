#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdint.h>
#include "gfparith.h"

void mpi_print(const char *c, const uint64_t *a, int len);
void mpi_conv_64to52(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_52to64(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_52to43(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_43to52(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_43to64(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_64to43(uint64_t *r, const uint64_t *a, int rlen, int alen);

__m512i set_vector(const uint64_t a7, const uint64_t a6, const uint64_t a5, const uint64_t a4, 
                   const uint64_t a3, const uint64_t a2, const uint64_t a1, const uint64_t a0);
void get_channel_8x1w(uint64_t *r, const htfe_t a, const int ch);

#endif