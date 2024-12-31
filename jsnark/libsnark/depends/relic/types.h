#ifndef TYPES_H
#define TYPES_H
#include <stdint.h>
#include <stdio.h>
#include "include/relic.h"
#include "padd.h"

#define NBITS 255                           // k's length
// #define WBITS 12                            // window's length
// #define NWINS ((NBITS + WBITS - 1) / WBITS) // number of windows  ((NBITS+WBITS-1)/WBITS)
// #define BUCKETNUM ((1 << WBITS) - 1)

#define DBL_P_SIZE  (1<<17)
#define AVX_WAY 8

typedef uint64_t bn52_t[HT_NWORDS];   

typedef struct
{
  bn52_t x;  // projective x coordinate in radix 2^52
  bn52_t y;  // projective y coordinate in radix 2^52
  bn52_t z;  // projective z coordinate in radix 2^52
} point52;

typedef point52 point52_t[1]; 

typedef bn52_t bn52_g2_t[2];

typedef struct
{
  bn52_g2_t x;  // projective x coordinate in radix 2^52
  bn52_g2_t y;  // projective y coordinate in radix 2^52
  bn52_g2_t z;  // projective z coordinate in radix 2^52
} point52_g2;

typedef point52_g2 point52_g2_t[1];

typedef struct
{
    int id;
    point52_t B;
    ep_t *wait;
} Bucket;

typedef struct
{
    int id;
    ep_t *P0;
    ep_t *P1;
} buf1_st;

typedef struct
{
    int id;
    point52_t point;
} point_buf;

typedef struct
{
    int id;
    point52_t *op1;
    point52_t *op2;
    int flag;
} buf2_st;

typedef struct
{
    ep2_t g2;
    ep_t g1;
} pair_st;

typedef struct
{
    point52_g2_t g2;
    point52_t g1;
} pair52_st;

typedef pair_st pair[1];
typedef pair52_st pair52[1];

typedef struct
{
    int id;
    pair52 B;
    pair *wait;
} pair_Bucket;

typedef struct
{
    int id;
    pair *P0;
    pair *P1;
}  pair_buf1_st;

typedef struct
{
    int id;
    pair52 point;
} pair_buf;

typedef struct
{
    int id;
    pair52 *op1;
    pair52 *op2;
    int flag;
} pair_buf2_st;

#endif