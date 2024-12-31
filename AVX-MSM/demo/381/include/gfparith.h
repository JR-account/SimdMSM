#ifndef _GFPARIT_H
#define _GFPARIT_H

#include <stdint.h>
#include "intrin.h"

// -----------------------------------------------------------------------------
// radix-52 parameters and constants for High-Throughput (HT) implementations

#define HT_BRADIX 52                    // limb size 
#define HT_NWORDS 8                    // limb number
#define HT_BMASK  0xFFFFFFFFFFFFFULL    // 2^52 - 1
#define HT_MONTW  0x3fffcfffcfffd //0x1301F632E294D       // the constant for Montgomery Multiplication

typedef __m512i htfe_t[HT_NWORDS];      // HT field element

// the prime p of the field
// static const uint64_t ht_p[HT_NWORDS] = {
//   0x1B90533C6C87B, 0xF457ACA8351B8, 0xF0B4F25C2721B, 0x5507516730CC1,
//   0xDA7AAC6C567F3, 0xFBFCC69322C9C, 0x83AEDC88C425A, 0x5E3E4C4AB42D0,
//   0xF89BFFC8AB0D1, 0x0065B48E8F740, };
static const uint64_t ht_p[HT_NWORDS] = {
  0xeffffffffaaab, 0xfeb153ffffb9f, 0x6b0f6241eabff, 0x12bf6730d2a0f,
  0x764774b84f385, 0x1ba7b6434bacd, 0x1ea397fe69a4b, 0x1a011, 
};


// p * 2
// static const uint64_t ht_pmul2[HT_NWORDS] = {
//   0x3720A678D90F6, 0xE8AF59506A370, 0xE169E4B84E437, 0xAA0EA2CE61983,
//   0xB4F558D8ACFE6, 0xF7F98D2645939, 0x075DB911884B5, 0xBC7C9895685A1,
//   0xF137FF91561A2, 0x00CB691D1EE81, };
static const uint64_t ht_pmul2[HT_NWORDS] = {
  0xdffffffff5556, 0xfd62a7ffff73f, 0xd61ec483d57ff, 0x257ece61a541e,
  0xec8ee9709e70a, 0x374f6c869759a, 0x3d472ffcd3496, 0x34022, 
};
 
// R mod p = 2^416 mod p
static const uint64_t ht_montR[HT_NWORDS] = {
  0x6480ea8e9b9af, 0x65766c8fe444f, 0x8b540fea96f7d, 0x3b2ee82efd422,
  0xa6723e5f0ade5, 0xff6eb6fdd4230, 0xe06ef23c24a25, 0x14c8e,
};  

// R^2 mod p = (2^520)^2 mod p 
// static const uint64_t ht_montR2[HT_NWORDS] = {
//   0x70C9A15C8CEBF, 0x5F05D4936EAAF, 0xF7EA4ADD4639B, 0x746FDEA7066EB,
//   0xE3AA0C3B19E9E, 0x3A0C5CE31A621, 0x114C5DF09803F, 0x702F11A883DAD,
//   0x94ECD5F5EC3F3, 0x00034FA8BE69F, };
// static const uint64_t ht_montR2[HT_NWORDS] = {
//   0x73281c647dac2, 0x94ace5aeaa152,
//   0x100af29d9c17b, 0xde4118a590a4c, 0x302d725a281d2, 0xc207df0fb68c2,
//   0xd270c3e5b7140, 0x8ca, 0x0, 0x0, };
// R^2 mod p = (2^416)^2 mod p 
static const uint64_t ht_montR2[HT_NWORDS] = {
  0xa5bf4cb89af51, 0x3afbba7ca31a2, 0x2646160ec71f1, 0xa84d710465903, 
  0x3480a4a188311, 0x98e5907ad91f5, 0x2075d74507266, 0x8746,  
};

// a = (2^384)^(-1) * (2^416)^2 mod p 
static const uint64_t mont_to_avx[HT_NWORDS] = {
  0x7fde37dba9366, 0x4e27525bc342b, 0x1f5b1e9778489, 0xb872b2b91b9dc, 
  0xb206f497dfcaf, 0x4137cc89a9b0b, 0xd9d20d7e39959, 0x411c,  
};

// static const uint64_t mont_to_avx[HT_NWORDS] = {
//   0x100000000, 0x0, 0x0, 0x0, 
//   0x0, 0x0, 0x0, 0x0,  
// };

// a =  (2^384) mod p 
static const uint64_t mont_to_relic[HT_NWORDS] = {
  0x900000002fffd, 0xbc40c0002760, 0x3c758baebf400, 0x57455f4898575, 
  0xd77ce58537052, 0x71a97a256ec6, 0xec3fa80e4935c, 0x15f65,  
};

// static const uint64_t mont_to_relic[HT_NWORDS] = {
//   0x28d21be8dc14d, 0x640583d131f06, 0x9587947f9ee17, 0x62ff5644a959, 
//   0x9e1dd6d0ff5, 0x8a18b41d687b8, 0x3b705b9a006ef, 0x61f,  
// };

// -----------------------------------------------------------------------------
// radix-64 constants 

// prime p in radix-64
static const uint64_t u64_p[8] = {
  0x1b81b90533c6c87B, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
  0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf, };

// (p-2) in radix-64
static const uint64_t u64_psub2[8] = {
  0x1b81b90533c6c879, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
  0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf, };

// (p-1)/2 in radix-64
static const uint64_t u64_pdiv2[8] = {
  0x8dc0dc8299e3643d, 0xe1390dfa2bd6541a, 0xa8b398660f85a792, 0xd3d56362b3f9aa83,
  0x2d7dfe63499164e6, 0x5a16841d76e44621, 0xfe455868af1f2625, 0x32da4747ba07c4df, };

// -----------------------------------------------------------------------------
// (8x1)-way prime-field operations
void gfp_add_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_dbl_8x1w(htfe_t r, const htfe_t a);
void gfp_sub_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_mul_8x1w(htfe_t r, const htfe_t a, const htfe_t b);
void gfp_sqr_8x1w(htfe_t r, const htfe_t a);
void gfp_rdcp_8x1w(htfe_t r, const htfe_t a);
void gfp_zero_8x1w(htfe_t r);
void gfp_mont_relic2avx_8x1w(htfe_t r, const htfe_t a);
void gfp_mont_avx2relic_8x1w(htfe_t r, const htfe_t a);
void gfp_num2mont_8x1w(htfe_t r, const htfe_t a);
void gfp_mont2num_8x1w(htfe_t r, const htfe_t a);
void gfp_copy_8x1w(htfe_t r, const htfe_t a);
int  gfp_iszero_8x1w(const htfe_t a);
#endif
