#ifndef PADD_H
#define PADD_H

#include "gfparith.h"
#include "gfp2arith.h"
#include "utils.h"

// the projective point on twisted Edwards curve with y and z coordinates 
typedef struct projective_point {
  htfe_t x;  // projective x coordinate
  htfe_t y;  // projective y coordinate
  htfe_t z;  // projective z coordinate
} htpoint;

typedef htpoint htpoint_t[1]; 

typedef struct projective_point2 {
  htfe2_t x;  // projective x coordinate
  htfe2_t y;  // projective y coordinate
  htfe2_t z;  // projective z coordinate
} ht2point;

typedef ht2point ht2point_t[1]; 

void PADD_jacob_8x1w_z1eqz2(htpoint_t r, const htpoint_t p, const htpoint_t q);
void PADD_jacob_8x1w_mix(htpoint_t r, const htpoint_t p, const htpoint_t q);
void PADD_jacob_8x1w_z1neqz2(htpoint_t r, const htpoint_t p, const htpoint_t q);

void PADD_projc_8x1w_z1eqz2(htpoint_t r, const htpoint_t p, const htpoint_t q);
void PADD_projc_8x1w_mix(htpoint_t r, const htpoint_t p, const htpoint_t q);
void PADD_projc_8x1w_z1neqz2(htpoint_t r, const htpoint_t p, const htpoint_t q);

void PADD_G2_projc_8x1w_z1eqz2(ht2point_t r, const ht2point_t p, const ht2point_t q);
void PADD_G2_projc_8x1w_mix(ht2point_t r, const ht2point_t p, const ht2point_t q);
void PADD_G2_projc_8x1w_z1neqz2(ht2point_t r, const ht2point_t p, const ht2point_t q);
#endif // PADD_H