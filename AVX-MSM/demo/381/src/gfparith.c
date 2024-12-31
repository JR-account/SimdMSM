#include "gfparith.h"

// (8x1)-way prime-field operations

// field addition r = a + b mod 2p
// a, b in [0, 2p) -> r in [0, 2p)
void gfp_add_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
	__m512i b5 = b[5], b6 = b[6], b7 = b[7];
	__m512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
	const __m512i vp0 = VSET1(ht_pmul2[0]), vp1 = VSET1(ht_pmul2[1]);
	const __m512i vp2 = VSET1(ht_pmul2[2]), vp3 = VSET1(ht_pmul2[3]);
	const __m512i vp4 = VSET1(ht_pmul2[4]), vp5 = VSET1(ht_pmul2[5]);
	const __m512i vp6 = VSET1(ht_pmul2[6]), vp7 = VSET1(ht_pmul2[7]);
	const __m512i vbmask = VSET1(HT_BMASK);

	// r = a + b
	r0 = VADD(a0, b0);
	r1 = VADD(a1, b1);
	r2 = VADD(a2, b2);
	r3 = VADD(a3, b3);
	r4 = VADD(a4, b4);
	r5 = VADD(a5, b5);
	r6 = VADD(a6, b6);
	r7 = VADD(a7, b7);

	// r = a + b - 2p
	r0 = VSUB(r0, vp0);
	r1 = VSUB(r1, vp1);
	r2 = VSUB(r2, vp2);
	r3 = VSUB(r3, vp3);
	r4 = VSUB(r4, vp4);
	r5 = VSUB(r5, vp5);
	r6 = VSUB(r6, vp6);
	r7 = VSUB(r7, vp7);

	// check the current r is positive or negative
	// carry propagation
	r1 = VADD(r1, VSRA(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSRA(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSRA(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSRA(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSRA(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSRA(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSRA(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);

	// if r is positive, then the corresponding element in smask = 0;
	// if r is negative, then the corresponding element in smask = all-1.
	smask = VSRA(r7, 63);
	// r = r + (2p & smask), add either 2p or 0 to the current r
	r0 = VADD(r0, VAND(vp0, smask));
	r1 = VADD(r1, VAND(vp1, smask));
	r2 = VADD(r2, VAND(vp2, smask));
	r3 = VADD(r3, VAND(vp3, smask));
	r4 = VADD(r4, VAND(vp4, smask));
	r5 = VADD(r5, VAND(vp5, smask));
	r6 = VADD(r6, VAND(vp6, smask));
	r7 = VADD(r7, VAND(vp7, smask));

	// carry propagation
	r1 = VADD(r1, VSHR(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSHR(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSHR(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSHR(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSHR(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSHR(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSHR(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);
	r7 = VAND(r7, vbmask);

	r[0] = r0;
	r[1] = r1;
	r[2] = r2;
	r[3] = r3;
	r[4] = r4;
	r[5] = r5;
	r[6] = r6;
	r[7] = r7;
}

// field addition r = a + b mod 2p
// a, b in [0, 2p) -> r in [0, 2p)
void gfp_dbl_8x1w(htfe_t r, const htfe_t a)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];

	__m512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
	const __m512i vp0 = VSET1(ht_pmul2[0]), vp1 = VSET1(ht_pmul2[1]);
	const __m512i vp2 = VSET1(ht_pmul2[2]), vp3 = VSET1(ht_pmul2[3]);
	const __m512i vp4 = VSET1(ht_pmul2[4]), vp5 = VSET1(ht_pmul2[5]);
	const __m512i vp6 = VSET1(ht_pmul2[6]), vp7 = VSET1(ht_pmul2[7]);
	const __m512i vbmask = VSET1(HT_BMASK);

	// r = a * 2
	r0 = VSHL(a0, 1);
	r1 = VSHL(a1, 1);
	r2 = VSHL(a2, 1);
	r3 = VSHL(a3, 1);
	r4 = VSHL(a4, 1);
	r5 = VSHL(a5, 1);
	r6 = VSHL(a6, 1);
	r7 = VSHL(a7, 1);

	// r = a * 2 - 2p
	r0 = VSUB(r0, vp0);
	r1 = VSUB(r1, vp1);
	r2 = VSUB(r2, vp2);
	r3 = VSUB(r3, vp3);
	r4 = VSUB(r4, vp4);
	r5 = VSUB(r5, vp5);
	r6 = VSUB(r6, vp6);
	r7 = VSUB(r7, vp7);

	// check the current r is positive or negative
	// carry propagation
	r1 = VADD(r1, VSRA(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSRA(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSRA(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSRA(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSRA(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSRA(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSRA(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);

	// if r is positive, then the corresponding element in smask = 0;
	// if r is negative, then the corresponding element in smask = all-1.
	smask = VSRA(r7, 63);
	// r = r + (2p & smask), add either 2p or 0 to the current r
	r0 = VADD(r0, VAND(vp0, smask));
	r1 = VADD(r1, VAND(vp1, smask));
	r2 = VADD(r2, VAND(vp2, smask));
	r3 = VADD(r3, VAND(vp3, smask));
	r4 = VADD(r4, VAND(vp4, smask));
	r5 = VADD(r5, VAND(vp5, smask));
	r6 = VADD(r6, VAND(vp6, smask));
	r7 = VADD(r7, VAND(vp7, smask));

	// carry propagation
	r1 = VADD(r1, VSHR(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSHR(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSHR(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSHR(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSHR(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSHR(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSHR(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);
	r7 = VAND(r7, vbmask);

	r[0] = r0;
	r[1] = r1;
	r[2] = r2;
	r[3] = r3;
	r[4] = r4;
	r[5] = r5;
	r[6] = r6;
	r[7] = r7;
}

// field subtraction r = a - b mod 2p
// a, b in [0, 2p) -> r in [0, 2p)
void gfp_sub_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
	__m512i b5 = b[5], b6 = b[6], b7 = b[7];
	__m512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
	const __m512i vp0 = VSET1(ht_pmul2[0]), vp1 = VSET1(ht_pmul2[1]);
	const __m512i vp2 = VSET1(ht_pmul2[2]), vp3 = VSET1(ht_pmul2[3]);
	const __m512i vp4 = VSET1(ht_pmul2[4]), vp5 = VSET1(ht_pmul2[5]);
	const __m512i vp6 = VSET1(ht_pmul2[6]), vp7 = VSET1(ht_pmul2[7]);
	const __m512i vbmask = VSET1(HT_BMASK);

	// r = a - b
	r0 = VSUB(a0, b0);
	r1 = VSUB(a1, b1);
	r2 = VSUB(a2, b2);
	r3 = VSUB(a3, b3);
	r4 = VSUB(a4, b4);
	r5 = VSUB(a5, b5);
	r6 = VSUB(a6, b6);
	r7 = VSUB(a7, b7);

	// check the current r is positive or negative
	// carry propagation
	r1 = VADD(r1, VSRA(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSRA(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSRA(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSRA(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSRA(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSRA(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSRA(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);

	// if r is positive, then the corresponding element in smask = 0;
	// if r is negative, then the corresponding element in smask = all-1.
	smask = VSRA(r7, 63);
	// r = r + (2p & smask), add either 2p or 0 to the current r
	r0 = VADD(r0, VAND(vp0, smask));
	r1 = VADD(r1, VAND(vp1, smask));
	r2 = VADD(r2, VAND(vp2, smask));
	r3 = VADD(r3, VAND(vp3, smask));
	r4 = VADD(r4, VAND(vp4, smask));
	r5 = VADD(r5, VAND(vp5, smask));
	r6 = VADD(r6, VAND(vp6, smask));
	r7 = VADD(r7, VAND(vp7, smask));

	// carry propagation
	r1 = VADD(r1, VSHR(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSHR(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSHR(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSHR(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSHR(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSHR(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSHR(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);
	r7 = VAND(r7, vbmask);

	r[0] = r0;
	r[1] = r1;
	r[2] = r2;
	r[3] = r3;
	r[4] = r4;
	r[5] = r5;
	r[6] = r6;
	r[7] = r7;
}

// Montgomery multiplication r = a * b mod 2p
// multiplication (product-scanning) interleaved with reduction (operand-scanning)
// -> r in [0, 2p)
void gfp_mul_8x1w(htfe_t r, const htfe_t a, const htfe_t b)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
	__m512i b5 = b[5], b6 = b[6], b7 = b[7];
	__m512i z0 = VZERO, z1 = VZERO, z2 = VZERO, z3 = VZERO, z4 = VZERO;
	__m512i z5 = VZERO, z6 = VZERO, z7 = VZERO, z8 = VZERO, z9 = VZERO;
	__m512i z10 = VZERO, z11 = VZERO, z12 = VZERO, z13 = VZERO, z14 = VZERO;
	__m512i z15 = VZERO;
	__m512i r0, r1, r2, r3, r4, r5, r6, r7, u;
	const __m512i vp0 = VSET1(ht_p[0]), vp1 = VSET1(ht_p[1]);
	const __m512i vp2 = VSET1(ht_p[2]), vp3 = VSET1(ht_p[3]);
	const __m512i vp4 = VSET1(ht_p[4]), vp5 = VSET1(ht_p[5]);
	const __m512i vp6 = VSET1(ht_p[6]), vp7 = VSET1(ht_p[7]);
	const __m512i vbmask = VSET1(HT_BMASK), vw = VSET1(HT_MONTW), zero = VZERO;

	// ---------------------------------------------------------------------------
	// 1st loop of integer multiplication

	z0 = VMACLO(z0, a0, b0);
	z1 = VMACHI(z1, a0, b0);

	z1 = VMACLO(z1, a0, b1);
	z1 = VMACLO(z1, a1, b0);
	z2 = VMACHI(z2, a0, b1);
	z2 = VMACHI(z2, a1, b0);

	z2 = VMACLO(z2, a0, b2);
	z2 = VMACLO(z2, a1, b1);
	z2 = VMACLO(z2, a2, b0);
	z3 = VMACHI(z3, a0, b2);
	z3 = VMACHI(z3, a1, b1);
	z3 = VMACHI(z3, a2, b0);

	z3 = VMACLO(z3, a0, b3);
	z3 = VMACLO(z3, a1, b2);
	z3 = VMACLO(z3, a2, b1);
	z3 = VMACLO(z3, a3, b0);
	z4 = VMACHI(z4, a0, b3);
	z4 = VMACHI(z4, a1, b2);
	z4 = VMACHI(z4, a2, b1);
	z4 = VMACHI(z4, a3, b0);

	z4 = VMACLO(z4, a0, b4);
	z4 = VMACLO(z4, a1, b3);
	z4 = VMACLO(z4, a2, b2);
	z4 = VMACLO(z4, a3, b1);
	z4 = VMACLO(z4, a4, b0);
	z5 = VMACHI(z5, a0, b4);
	z5 = VMACHI(z5, a1, b3);
	z5 = VMACHI(z5, a2, b2);
	z5 = VMACHI(z5, a3, b1);
	z5 = VMACHI(z5, a4, b0);

	z5 = VMACLO(z5, a0, b5);
	z5 = VMACLO(z5, a1, b4);
	z5 = VMACLO(z5, a2, b3);
	z5 = VMACLO(z5, a3, b2);
	z5 = VMACLO(z5, a4, b1);
	z5 = VMACLO(z5, a5, b0);
	z6 = VMACHI(z6, a0, b5);
	z6 = VMACHI(z6, a1, b4);
	z6 = VMACHI(z6, a2, b3);
	z6 = VMACHI(z6, a3, b2);
	z6 = VMACHI(z6, a4, b1);
	z6 = VMACHI(z6, a5, b0);

	z6 = VMACLO(z6, a0, b6);
	z6 = VMACLO(z6, a1, b5);
	z6 = VMACLO(z6, a2, b4);
	z6 = VMACLO(z6, a3, b3);
	z6 = VMACLO(z6, a4, b2);
	z6 = VMACLO(z6, a5, b1);
	z6 = VMACLO(z6, a6, b0);
	z7 = VMACHI(z7, a0, b6);
	z7 = VMACHI(z7, a1, b5);
	z7 = VMACHI(z7, a2, b4);
	z7 = VMACHI(z7, a3, b3);
	z7 = VMACHI(z7, a4, b2);
	z7 = VMACHI(z7, a5, b1);
	z7 = VMACHI(z7, a6, b0);

	z7 = VMACLO(z7, a0, b7);
	z7 = VMACLO(z7, a1, b6);
	z7 = VMACLO(z7, a2, b5);
	z7 = VMACLO(z7, a3, b4);
	z7 = VMACLO(z7, a4, b3);
	z7 = VMACLO(z7, a5, b2);
	z7 = VMACLO(z7, a6, b1);
	z7 = VMACLO(z7, a7, b0);
	z8 = VMACHI(z8, a0, b7);
	z8 = VMACHI(z8, a1, b6);
	z8 = VMACHI(z8, a2, b5);
	z8 = VMACHI(z8, a3, b4);
	z8 = VMACHI(z8, a4, b3);
	z8 = VMACHI(z8, a5, b2);
	z8 = VMACHI(z8, a6, b1);
	z8 = VMACHI(z8, a7, b0);

	// ---------------------------------------------------------------------------
	// 2nd loop of integer multiplication + Montgomery reduction

	u = VMACLO(zero, z0, vw);
	z0 = VMACLO(z0, u, vp0);
	z1 = VMACHI(z1, u, vp0);
	z1 = VMACLO(z1, u, vp1);
	z2 = VMACHI(z2, u, vp1);
	z2 = VMACLO(z2, u, vp2);
	z3 = VMACHI(z3, u, vp2);
	z3 = VMACLO(z3, u, vp3);
	z4 = VMACHI(z4, u, vp3);
	z4 = VMACLO(z4, u, vp4);
	z5 = VMACHI(z5, u, vp4);
	z5 = VMACLO(z5, u, vp5);
	z6 = VMACHI(z6, u, vp5);
	z6 = VMACLO(z6, u, vp6);
	z7 = VMACHI(z7, u, vp6);
	z7 = VMACLO(z7, u, vp7);
	z8 = VMACHI(z8, u, vp7);
	z1 = VADD(z1, VSHR(z0, HT_BRADIX));

	z8 = VMACLO(z8, a1, b7);
	z8 = VMACLO(z8, a2, b6);
	z8 = VMACLO(z8, a3, b5);
	z8 = VMACLO(z8, a4, b4);
	z8 = VMACLO(z8, a5, b3);
	z8 = VMACLO(z8, a6, b2);
	z8 = VMACLO(z8, a7, b1);
	z9 = VMACHI(z9, a1, b7);
	z9 = VMACHI(z9, a2, b6);
	z9 = VMACHI(z9, a3, b5);
	z9 = VMACHI(z9, a4, b4);
	z9 = VMACHI(z9, a5, b3);
	z9 = VMACHI(z9, a6, b2);
	z9 = VMACHI(z9, a7, b1);

	u = VMACLO(zero, z1, vw);
	z1 = VMACLO(z1, u, vp0);
	z2 = VMACHI(z2, u, vp0);
	z2 = VMACLO(z2, u, vp1);
	z3 = VMACHI(z3, u, vp1);
	z3 = VMACLO(z3, u, vp2);
	z4 = VMACHI(z4, u, vp2);
	z4 = VMACLO(z4, u, vp3);
	z5 = VMACHI(z5, u, vp3);
	z5 = VMACLO(z5, u, vp4);
	z6 = VMACHI(z6, u, vp4);
	z6 = VMACLO(z6, u, vp5);
	z7 = VMACHI(z7, u, vp5);
	z7 = VMACLO(z7, u, vp6);
	z8 = VMACHI(z8, u, vp6);
	z8 = VMACLO(z8, u, vp7);
	z9 = VMACHI(z9, u, vp7);
	z2 = VADD(z2, VSHR(z1, HT_BRADIX));

	z9 = VMACLO(z9, a2, b7);
	z9 = VMACLO(z9, a3, b6);
	z9 = VMACLO(z9, a4, b5);
	z9 = VMACLO(z9, a5, b4);
	z9 = VMACLO(z9, a6, b3);
	z9 = VMACLO(z9, a7, b2);
	z10 = VMACHI(z10, a2, b7);
	z10 = VMACHI(z10, a3, b6);
	z10 = VMACHI(z10, a4, b5);
	z10 = VMACHI(z10, a5, b4);
	z10 = VMACHI(z10, a6, b3);
	z10 = VMACHI(z10, a7, b2);

	u = VMACLO(zero, z2, vw);
	z2 = VMACLO(z2, u, vp0);
	z3 = VMACHI(z3, u, vp0);
	z3 = VMACLO(z3, u, vp1);
	z4 = VMACHI(z4, u, vp1);
	z4 = VMACLO(z4, u, vp2);
	z5 = VMACHI(z5, u, vp2);
	z5 = VMACLO(z5, u, vp3);
	z6 = VMACHI(z6, u, vp3);
	z6 = VMACLO(z6, u, vp4);
	z7 = VMACHI(z7, u, vp4);
	z7 = VMACLO(z7, u, vp5);
	z8 = VMACHI(z8, u, vp5);
	z8 = VMACLO(z8, u, vp6);
	z9 = VMACHI(z9, u, vp6);
	z9 = VMACLO(z9, u, vp7);
	z10 = VMACHI(z10, u, vp7);
	z3 = VADD(z3, VSHR(z2, HT_BRADIX));

	z10 = VMACLO(z10, a3, b7);
	z10 = VMACLO(z10, a4, b6);
	z10 = VMACLO(z10, a5, b5);
	z10 = VMACLO(z10, a6, b4);
	z10 = VMACLO(z10, a7, b3);
	z11 = VMACHI(z11, a3, b7);
	z11 = VMACHI(z11, a4, b6);
	z11 = VMACHI(z11, a5, b5);
	z11 = VMACHI(z11, a6, b4);
	z11 = VMACHI(z11, a7, b3);

	u = VMACLO(zero, z3, vw);
	z3 = VMACLO(z3, u, vp0);
	z4 = VMACHI(z4, u, vp0);
	z4 = VMACLO(z4, u, vp1);
	z5 = VMACHI(z5, u, vp1);
	z5 = VMACLO(z5, u, vp2);
	z6 = VMACHI(z6, u, vp2);
	z6 = VMACLO(z6, u, vp3);
	z7 = VMACHI(z7, u, vp3);
	z7 = VMACLO(z7, u, vp4);
	z8 = VMACHI(z8, u, vp4);
	z8 = VMACLO(z8, u, vp5);
	z9 = VMACHI(z9, u, vp5);
	z9 = VMACLO(z9, u, vp6);
	z10 = VMACHI(z10, u, vp6);
	z10 = VMACLO(z10, u, vp7);
	z11 = VMACHI(z11, u, vp7);
	z4 = VADD(z4, VSHR(z3, HT_BRADIX));

	z11 = VMACLO(z11, a4, b7);
	z11 = VMACLO(z11, a5, b6);
	z11 = VMACLO(z11, a6, b5);
	z11 = VMACLO(z11, a7, b4);
	z12 = VMACHI(z12, a4, b7);
	z12 = VMACHI(z12, a5, b6);
	z12 = VMACHI(z12, a6, b5);
	z12 = VMACHI(z12, a7, b4);

	u = VMACLO(zero, z4, vw);
	z4 = VMACLO(z4, u, vp0);
	z5 = VMACHI(z5, u, vp0);
	z5 = VMACLO(z5, u, vp1);
	z6 = VMACHI(z6, u, vp1);
	z6 = VMACLO(z6, u, vp2);
	z7 = VMACHI(z7, u, vp2);
	z7 = VMACLO(z7, u, vp3);
	z8 = VMACHI(z8, u, vp3);
	z8 = VMACLO(z8, u, vp4);
	z9 = VMACHI(z9, u, vp4);
	z9 = VMACLO(z9, u, vp5);
	z10 = VMACHI(z10, u, vp5);
	z10 = VMACLO(z10, u, vp6);
	z11 = VMACHI(z11, u, vp6);
	z11 = VMACLO(z11, u, vp7);
	z12 = VMACHI(z12, u, vp7);
	z5 = VADD(z5, VSHR(z4, HT_BRADIX));

	z12 = VMACLO(z12, a5, b7);
	z12 = VMACLO(z12, a6, b6);
	z12 = VMACLO(z12, a7, b5);
	z13 = VMACHI(z13, a5, b7);
	z13 = VMACHI(z13, a6, b6);
	z13 = VMACHI(z13, a7, b5);

	u = VMACLO(zero, z5, vw);
	z5 = VMACLO(z5, u, vp0);
	z6 = VMACHI(z6, u, vp0);
	z6 = VMACLO(z6, u, vp1);
	z7 = VMACHI(z7, u, vp1);
	z7 = VMACLO(z7, u, vp2);
	z8 = VMACHI(z8, u, vp2);
	z8 = VMACLO(z8, u, vp3);
	z9 = VMACHI(z9, u, vp3);
	z9 = VMACLO(z9, u, vp4);
	z10 = VMACHI(z10, u, vp4);
	z10 = VMACLO(z10, u, vp5);
	z11 = VMACHI(z11, u, vp5);
	z11 = VMACLO(z11, u, vp6);
	z12 = VMACHI(z12, u, vp6);
	z12 = VMACLO(z12, u, vp7);
	z13 = VMACHI(z13, u, vp7);
	z6 = VADD(z6, VSHR(z5, HT_BRADIX));

	z13 = VMACLO(z13, a6, b7);
	z13 = VMACLO(z13, a7, b6);
	z14 = VMACHI(z14, a6, b7);
	z14 = VMACHI(z14, a7, b6);

	u = VMACLO(zero, z6, vw);
	z6 = VMACLO(z6, u, vp0);
	z7 = VMACHI(z7, u, vp0);
	z7 = VMACLO(z7, u, vp1);
	z8 = VMACHI(z8, u, vp1);
	z8 = VMACLO(z8, u, vp2);
	z9 = VMACHI(z9, u, vp2);
	z9 = VMACLO(z9, u, vp3);
	z10 = VMACHI(z10, u, vp3);
	z10 = VMACLO(z10, u, vp4);
	z11 = VMACHI(z11, u, vp4);
	z11 = VMACLO(z11, u, vp5);
	z12 = VMACHI(z12, u, vp5);
	z12 = VMACLO(z12, u, vp6);
	z13 = VMACHI(z13, u, vp6);
	z13 = VMACLO(z13, u, vp7);
	z14 = VMACHI(z14, u, vp7);
	z7 = VADD(z7, VSHR(z6, HT_BRADIX));

	z14 = VMACLO(z14, a7, b7);
	z15 = VMACHI(z15, a7, b7);

	u = VMACLO(zero, z7, vw);
	z7 = VMACLO(z7, u, vp0);
	z8 = VMACHI(z8, u, vp0);
	z8 = VMACLO(z8, u, vp1);
	z9 = VMACHI(z9, u, vp1);
	z9 = VMACLO(z9, u, vp2);
	z10 = VMACHI(z10, u, vp2);
	z10 = VMACLO(z10, u, vp3);
	z11 = VMACHI(z11, u, vp3);
	z11 = VMACLO(z11, u, vp4);
	z12 = VMACHI(z12, u, vp4);
	z12 = VMACLO(z12, u, vp5);
	z13 = VMACHI(z13, u, vp5);
	z13 = VMACLO(z13, u, vp6);
	z14 = VMACHI(z14, u, vp6);
	z14 = VMACLO(z14, u, vp7);
	z15 = VMACHI(z15, u, vp7);
	z8 = VADD(z8, VSHR(z7, HT_BRADIX));

	// carry propagation
	z9 = VADD(z9, VSHR(z8, HT_BRADIX));
	z8 = VAND(z8, vbmask);
	z10 = VADD(z10, VSHR(z9, HT_BRADIX));
	z9 = VAND(z9, vbmask);
	z11 = VADD(z11, VSHR(z10, HT_BRADIX));
	z10 = VAND(z10, vbmask);
	z12 = VADD(z12, VSHR(z11, HT_BRADIX));
	z11 = VAND(z11, vbmask);
	z13 = VADD(z13, VSHR(z12, HT_BRADIX));
	z12 = VAND(z12, vbmask);
	z14 = VADD(z14, VSHR(z13, HT_BRADIX));
	z13 = VAND(z13, vbmask);
	z15 = VADD(z15, VSHR(z14, HT_BRADIX));
	z14 = VAND(z14, vbmask);

	// ---------------------------------------------------------------------------

	r[0] = z8;
	r[1] = z9;
	r[2] = z10;
	r[3] = z11;
	r[4] = z12;
	r[5] = z13;
	r[6] = z14;
	r[7] = z15;
}

// Montgomery squaring r = a^2 mod 2p
// squaring (product-scanning) interleaved with reduction (operand-scanning)
// -> r in [0, 2p)
void gfp_sqr_8x1w(htfe_t r, const htfe_t a)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i z0 = VZERO, z1 = VZERO, z2 = VZERO, z3 = VZERO, z4 = VZERO;
	__m512i z5 = VZERO, z6 = VZERO, z7 = VZERO, z8 = VZERO, z9 = VZERO;
	__m512i z10 = VZERO, z11 = VZERO, z12 = VZERO, z13 = VZERO, z14 = VZERO;
	__m512i z15 = VZERO;
	__m512i r0, r1, r2, r3, r4, r5, r6, r7, u;
	const __m512i vp0 = VSET1(ht_p[0]), vp1 = VSET1(ht_p[1]);
	const __m512i vp2 = VSET1(ht_p[2]), vp3 = VSET1(ht_p[3]);
	const __m512i vp4 = VSET1(ht_p[4]), vp5 = VSET1(ht_p[5]);
	const __m512i vp6 = VSET1(ht_p[6]), vp7 = VSET1(ht_p[7]);
	const __m512i vbmask = VSET1(HT_BMASK), vw = VSET1(HT_MONTW), zero = VZERO;

	// ---------------------------------------------------------------------------
	// 1st loop of integer squaring

	z1 = VMACLO(z1, a0, a1);
	z2 = VMACHI(z2, a0, a1);
	z1 = VADD(z1, z1);
	z0 = VMACLO(z0, a0, a0);
	z1 = VMACHI(z1, a0, a0);

	z2 = VMACLO(z2, a0, a2);
	z3 = VMACHI(z3, a0, a2);
	z2 = VADD(z2, z2);

	z3 = VMACLO(z3, a0, a3);
	z3 = VMACLO(z3, a1, a2);
	z4 = VMACHI(z4, a0, a3);
	z4 = VMACHI(z4, a1, a2);
	z3 = VADD(z3, z3);
	z2 = VMACLO(z2, a1, a1);
	z3 = VMACHI(z3, a1, a1);

	z4 = VMACLO(z4, a0, a4);
	z4 = VMACLO(z4, a1, a3);
	z5 = VMACHI(z5, a0, a4);
	z5 = VMACHI(z5, a1, a3);
	z4 = VADD(z4, z4);

	z5 = VMACLO(z5, a0, a5);
	z5 = VMACLO(z5, a1, a4);
	z5 = VMACLO(z5, a2, a3);
	z6 = VMACHI(z6, a0, a5);
	z6 = VMACHI(z6, a1, a4);
	z6 = VMACHI(z6, a2, a3);
	z5 = VADD(z5, z5);
	z4 = VMACLO(z4, a2, a2);
	z5 = VMACHI(z5, a2, a2);

	z6 = VMACLO(z6, a0, a6);
	z6 = VMACLO(z6, a1, a5);
	z6 = VMACLO(z6, a2, a4);
	z7 = VMACHI(z7, a0, a6);
	z7 = VMACHI(z7, a1, a5);
	z7 = VMACHI(z7, a2, a4);
	z6 = VADD(z6, z6);

	z7 = VMACLO(z7, a0, a7);
	z7 = VMACLO(z7, a1, a6);
	z7 = VMACLO(z7, a2, a5);
	z7 = VMACLO(z7, a3, a4);
	z8 = VMACHI(z8, a0, a7);
	z8 = VMACHI(z8, a1, a6);
	z8 = VMACHI(z8, a2, a5);
	z8 = VMACHI(z8, a3, a4);
	z7 = VADD(z7, z7);
	z6 = VMACLO(z6, a3, a3);
	z7 = VMACHI(z7, a3, a3);

	// ---------------------------------------------------------------------------
	// 2nd loop of integer squaring + Montgomery reduction

	z8 = VMACLO(z8, a1, a7);
	z8 = VMACLO(z8, a2, a6);
	z8 = VMACLO(z8, a3, a5);
	z9 = VMACHI(z9, a1, a7);
	z9 = VMACHI(z9, a2, a6);
	z9 = VMACHI(z9, a3, a5);
	z8 = VADD(z8, z8);

	u = VMACLO(zero, z0, vw);
	z0 = VMACLO(z0, u, vp0);
	z1 = VMACHI(z1, u, vp0);
	z1 = VMACLO(z1, u, vp1);
	z2 = VMACHI(z2, u, vp1);
	z2 = VMACLO(z2, u, vp2);
	z3 = VMACHI(z3, u, vp2);
	z3 = VMACLO(z3, u, vp3);
	z4 = VMACHI(z4, u, vp3);
	z4 = VMACLO(z4, u, vp4);
	z5 = VMACHI(z5, u, vp4);
	z5 = VMACLO(z5, u, vp5);
	z6 = VMACHI(z6, u, vp5);
	z6 = VMACLO(z6, u, vp6);
	z7 = VMACHI(z7, u, vp6);
	z7 = VMACLO(z7, u, vp7);
	z8 = VMACHI(z8, u, vp7);
	z1 = VADD(z1, VSHR(z0, HT_BRADIX));

	z9 = VMACLO(z9, a2, a7);
	z9 = VMACLO(z9, a3, a6);
	z9 = VMACLO(z9, a4, a5);
	z10 = VMACHI(z10, a2, a7);
	z10 = VMACHI(z10, a3, a6);
	z10 = VMACHI(z10, a4, a5);
	z9 = VADD(z9, z9);
	z8 = VMACLO(z8, a4, a4);
	z9 = VMACHI(z9, a4, a4);

	u = VMACLO(zero, z1, vw);
	z1 = VMACLO(z1, u, vp0);
	z2 = VMACHI(z2, u, vp0);
	z2 = VMACLO(z2, u, vp1);
	z3 = VMACHI(z3, u, vp1);
	z3 = VMACLO(z3, u, vp2);
	z4 = VMACHI(z4, u, vp2);
	z4 = VMACLO(z4, u, vp3);
	z5 = VMACHI(z5, u, vp3);
	z5 = VMACLO(z5, u, vp4);
	z6 = VMACHI(z6, u, vp4);
	z6 = VMACLO(z6, u, vp5);
	z7 = VMACHI(z7, u, vp5);
	z7 = VMACLO(z7, u, vp6);
	z8 = VMACHI(z8, u, vp6);
	z8 = VMACLO(z8, u, vp7);
	z9 = VMACHI(z9, u, vp7);
	z2 = VADD(z2, VSHR(z1, HT_BRADIX));

	z10 = VMACLO(z10, a3, a7);
	z10 = VMACLO(z10, a4, a6);
	z11 = VMACHI(z11, a3, a7);
	z11 = VMACHI(z11, a4, a6);
	z10 = VADD(z10, z10);

	u = VMACLO(zero, z2, vw);
	z2 = VMACLO(z2, u, vp0);
	z3 = VMACHI(z3, u, vp0);
	z3 = VMACLO(z3, u, vp1);
	z4 = VMACHI(z4, u, vp1);
	z4 = VMACLO(z4, u, vp2);
	z5 = VMACHI(z5, u, vp2);
	z5 = VMACLO(z5, u, vp3);
	z6 = VMACHI(z6, u, vp3);
	z6 = VMACLO(z6, u, vp4);
	z7 = VMACHI(z7, u, vp4);
	z7 = VMACLO(z7, u, vp5);
	z8 = VMACHI(z8, u, vp5);
	z8 = VMACLO(z8, u, vp6);
	z9 = VMACHI(z9, u, vp6);
	z9 = VMACLO(z9, u, vp7);
	z10 = VMACHI(z10, u, vp7);
	z3 = VADD(z3, VSHR(z2, HT_BRADIX));

	z11 = VMACLO(z11, a4, a7);
	z11 = VMACLO(z11, a5, a6);
	z12 = VMACHI(z12, a4, a7);
	z12 = VMACHI(z12, a5, a6);
	z11 = VADD(z11, z11);
	z10 = VMACLO(z10, a5, a5);
	z11 = VMACHI(z11, a5, a5);

	u = VMACLO(zero, z3, vw);
	z3 = VMACLO(z3, u, vp0);
	z4 = VMACHI(z4, u, vp0);
	z4 = VMACLO(z4, u, vp1);
	z5 = VMACHI(z5, u, vp1);
	z5 = VMACLO(z5, u, vp2);
	z6 = VMACHI(z6, u, vp2);
	z6 = VMACLO(z6, u, vp3);
	z7 = VMACHI(z7, u, vp3);
	z7 = VMACLO(z7, u, vp4);
	z8 = VMACHI(z8, u, vp4);
	z8 = VMACLO(z8, u, vp5);
	z9 = VMACHI(z9, u, vp5);
	z9 = VMACLO(z9, u, vp6);
	z10 = VMACHI(z10, u, vp6);
	z10 = VMACLO(z10, u, vp7);
	z11 = VMACHI(z11, u, vp7);
	z4 = VADD(z4, VSHR(z3, HT_BRADIX));

	z12 = VMACLO(z12, a5, a7);
	z13 = VMACHI(z13, a5, a7);
	z12 = VADD(z12, z12);

	u = VMACLO(zero, z4, vw);
	z4 = VMACLO(z4, u, vp0);
	z5 = VMACHI(z5, u, vp0);
	z5 = VMACLO(z5, u, vp1);
	z6 = VMACHI(z6, u, vp1);
	z6 = VMACLO(z6, u, vp2);
	z7 = VMACHI(z7, u, vp2);
	z7 = VMACLO(z7, u, vp3);
	z8 = VMACHI(z8, u, vp3);
	z8 = VMACLO(z8, u, vp4);
	z9 = VMACHI(z9, u, vp4);
	z9 = VMACLO(z9, u, vp5);
	z10 = VMACHI(z10, u, vp5);
	z10 = VMACLO(z10, u, vp6);
	z11 = VMACHI(z11, u, vp6);
	z11 = VMACLO(z11, u, vp7);
	z12 = VMACHI(z12, u, vp7);
	z5 = VADD(z5, VSHR(z4, HT_BRADIX));

	z13 = VMACLO(z13, a6, a7);
	z14 = VMACHI(z14, a6, a7);
	z13 = VADD(z13, z13);
	z12 = VMACLO(z12, a6, a6);
	z13 = VMACHI(z13, a6, a6);

	u = VMACLO(zero, z5, vw);
	z5 = VMACLO(z5, u, vp0);
	z6 = VMACHI(z6, u, vp0);
	z6 = VMACLO(z6, u, vp1);
	z7 = VMACHI(z7, u, vp1);
	z7 = VMACLO(z7, u, vp2);
	z8 = VMACHI(z8, u, vp2);
	z8 = VMACLO(z8, u, vp3);
	z9 = VMACHI(z9, u, vp3);
	z9 = VMACLO(z9, u, vp4);
	z10 = VMACHI(z10, u, vp4);
	z10 = VMACLO(z10, u, vp5);
	z11 = VMACHI(z11, u, vp5);
	z11 = VMACLO(z11, u, vp6);
	z12 = VMACHI(z12, u, vp6);
	z12 = VMACLO(z12, u, vp7);
	z13 = VMACHI(z13, u, vp7);
	z6 = VADD(z6, VSHR(z5, HT_BRADIX));

	z14 = VADD(z14, z14);

	u = VMACLO(zero, z6, vw);
	z6 = VMACLO(z6, u, vp0);
	z7 = VMACHI(z7, u, vp0);
	z7 = VMACLO(z7, u, vp1);
	z8 = VMACHI(z8, u, vp1);
	z8 = VMACLO(z8, u, vp2);
	z9 = VMACHI(z9, u, vp2);
	z9 = VMACLO(z9, u, vp3);
	z10 = VMACHI(z10, u, vp3);
	z10 = VMACLO(z10, u, vp4);
	z11 = VMACHI(z11, u, vp4);
	z11 = VMACLO(z11, u, vp5);
	z12 = VMACHI(z12, u, vp5);
	z12 = VMACLO(z12, u, vp6);
	z13 = VMACHI(z13, u, vp6);
	z13 = VMACLO(z13, u, vp7);
	z14 = VMACHI(z14, u, vp7);
	z7 = VADD(z7, VSHR(z6, HT_BRADIX));

	z14 = VMACLO(z14, a7, a7);
	z15 = VMACHI(z15, a7, a7);

	u = VMACLO(zero, z7, vw);
	z7 = VMACLO(z7, u, vp0);
	z8 = VMACHI(z8, u, vp0);
	z8 = VMACLO(z8, u, vp1);
	z9 = VMACHI(z9, u, vp1);
	z9 = VMACLO(z9, u, vp2);
	z10 = VMACHI(z10, u, vp2);
	z10 = VMACLO(z10, u, vp3);
	z11 = VMACHI(z11, u, vp3);
	z11 = VMACLO(z11, u, vp4);
	z12 = VMACHI(z12, u, vp4);
	z12 = VMACLO(z12, u, vp5);
	z13 = VMACHI(z13, u, vp5);
	z13 = VMACLO(z13, u, vp6);
	z14 = VMACHI(z14, u, vp6);
	z14 = VMACLO(z14, u, vp7);
	z15 = VMACHI(z15, u, vp7);
	z8 = VADD(z8, VSHR(z7, HT_BRADIX));

	z9 = VADD(z9, VSHR(z8, HT_BRADIX));
	z8 = VAND(z8, vbmask);
	z10 = VADD(z10, VSHR(z9, HT_BRADIX));
	z9 = VAND(z9, vbmask);
	z11 = VADD(z11, VSHR(z10, HT_BRADIX));
	z10 = VAND(z10, vbmask);
	z12 = VADD(z12, VSHR(z11, HT_BRADIX));
	z11 = VAND(z11, vbmask);
	z13 = VADD(z13, VSHR(z12, HT_BRADIX));
	z12 = VAND(z12, vbmask);
	z14 = VADD(z14, VSHR(z13, HT_BRADIX));
	z13 = VAND(z13, vbmask);
	z15 = VADD(z15, VSHR(z14, HT_BRADIX));
	z14 = VAND(z14, vbmask);

	// ---------------------------------------------------------------------------

	r[0] = z8;
	r[1] = z9;
	r[2] = z10;
	r[3] = z11;
	r[4] = z12;
	r[5] = z13;
	r[6] = z14;
	r[7] = z15;
}

// reduce the field element from [0, 2p) to [0, p)
// a in [0, 2p) -> r in [0, p)
void gfp_rdcp_8x1w(htfe_t r, const htfe_t a)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
	const __m512i vp0 = VSET1(ht_p[0]), vp1 = VSET1(ht_p[1]);
	const __m512i vp2 = VSET1(ht_p[2]), vp3 = VSET1(ht_p[3]);
	const __m512i vp4 = VSET1(ht_p[4]), vp5 = VSET1(ht_p[5]);
	const __m512i vp6 = VSET1(ht_p[6]), vp7 = VSET1(ht_p[7]);
	const __m512i vbmask = VSET1(HT_BMASK);

	// r = a - p
	r0 = VSUB(a0, vp0);
	r1 = VSUB(a1, vp1);
	r2 = VSUB(a2, vp2);
	r3 = VSUB(a3, vp3);
	r4 = VSUB(a4, vp4);
	r5 = VSUB(a5, vp5);
	r6 = VSUB(a6, vp6);
	r7 = VSUB(a7, vp7);

	// check the current r is positive or negative
	// carry propagation
	r1 = VADD(r1, VSRA(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSRA(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSRA(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSRA(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSRA(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSRA(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSRA(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);

	// if r is positive, then the corresponding element in smask = 0;
	// if r is negative, then the corresponding element in smask = all-1.
	smask = VSRA(r7, 63);
	// r = r + (p & smask), add either p or 0 to the current r
	r0 = VADD(r0, VAND(vp0, smask));
	r1 = VADD(r1, VAND(vp1, smask));
	r2 = VADD(r2, VAND(vp2, smask));
	r3 = VADD(r3, VAND(vp3, smask));
	r4 = VADD(r4, VAND(vp4, smask));
	r5 = VADD(r5, VAND(vp5, smask));
	r6 = VADD(r6, VAND(vp6, smask));
	r7 = VADD(r7, VAND(vp7, smask));

	// carry propagation
	r1 = VADD(r1, VSHR(r0, HT_BRADIX));
	r0 = VAND(r0, vbmask);
	r2 = VADD(r2, VSHR(r1, HT_BRADIX));
	r1 = VAND(r1, vbmask);
	r3 = VADD(r3, VSHR(r2, HT_BRADIX));
	r2 = VAND(r2, vbmask);
	r4 = VADD(r4, VSHR(r3, HT_BRADIX));
	r3 = VAND(r3, vbmask);
	r5 = VADD(r5, VSHR(r4, HT_BRADIX));
	r4 = VAND(r4, vbmask);
	r6 = VADD(r6, VSHR(r5, HT_BRADIX));
	r5 = VAND(r5, vbmask);
	r7 = VADD(r7, VSHR(r6, HT_BRADIX));
	r6 = VAND(r6, vbmask);
	r7 = VAND(r7, vbmask);

	r[0] = r0;
	r[1] = r1;
	r[2] = r2;
	r[3] = r3;
	r[4] = r4;
	r[5] = r5;
	r[6] = r6;
	r[7] = r7;
}

// set the field element r to be 0
// -> r = 0
void gfp_zero_8x1w(htfe_t r)
{
	r[0] = VZERO;
	r[1] = VZERO;
	r[2] = VZERO;
	r[3] = VZERO;
	r[4] = VZERO;
	r[5] = VZERO;
	r[6] = VZERO;
	r[7] = VZERO;
}

// convert an integer from number domain to Montgomery domain r = a * R^2 mod 2p
// -> r in [0, 2p)
void gfp_num2mont_8x1w(htfe_t r, const htfe_t a)
{
	int i;
	htfe_t vR2;

	for (i = 0; i < HT_NWORDS; i++)
		vR2[i] = VSET1(ht_montR2[i]);
	gfp_mul_8x1w(r, a, vR2); // r = a * R^2 mod 2p
}

// convert an integer from Montgomery domain 2^384 to Montgomery domain 2^416 r = aR1 * R1^(-1)R2^2 mod 2p
// -> r in [0, 2p)
void gfp_mont_relic2avx_8x1w(htfe_t r, const htfe_t a)
{
	int i;
	htfe_t vR2;

	for (i = 0; i < HT_NWORDS; i++)
		vR2[i] = VSET1(mont_to_avx[i]);
	gfp_mul_8x1w(r, a, vR2);
}

// convert an integer from Montgomery domain 2^416 to Montgomery domain 2^384 r = aR1 * R1^(-1)R2^2 mod 2p
// -> r in [0, 2p)
void gfp_mont_avx2relic_8x1w(htfe_t r, const htfe_t a)
{
	int i;
	htfe_t vR2;

	for (i = 0; i < HT_NWORDS; i++)
		vR2[i] = VSET1(mont_to_relic[i]);
	gfp_mul_8x1w(r, a, vR2);
	// gfp_rdcp_8x1w(r, r);	  // r = r mod p
}

// convert an integer from Montgomery domain to number domain r = a * 1 mod 2p
// -> r in [0, p)
void gfp_mont2num_8x1w(htfe_t r, const htfe_t a)
{
	htfe_t vone;

	gfp_zero_8x1w(vone);	  // vone = 0
	vone[0] = VSET1(1);		  // vone = 1
	gfp_mul_8x1w(r, a, vone); // r = a * 1 mod 2p
	gfp_rdcp_8x1w(r, r);	  // r = r mod p
}

// copy field element a to r
void gfp_copy_8x1w(htfe_t r, const htfe_t a)
{
	r[0] = a[0];
	r[1] = a[1];
	r[2] = a[2];
	r[3] = a[3];
	r[4] = a[4];
	r[5] = a[5];
	r[6] = a[6];
	r[7] = a[7];
}

// check whether the field element is zero
// if a is     zero, return 1;
// if a is not zero, return 0.
int gfp_iszero_8x1w(const htfe_t a)
{
	__m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	__m512i a5 = a[5], a6 = a[6], a7 = a[7];
	__m512i t = VZERO;

	// t = OR all the limbs of a
	t = VOR(t, a0);
	t = VOR(t, a1);
	t = VOR(t, a2);
	t = VOR(t, a3);
	t = VOR(t, a4);
	t = VOR(t, a5);
	t = VOR(t, a6);
	t = VOR(t, a7);

	// if a == 0 then t is all-0; otherwise t is non-0
	// t = VAND(t, VSET1(HT_BMASK)); // t is 52-bit now
	// t = VSUB(VZERO, t);			  // t is all-0 (if a is zero) or negative (if a is non-zero)
	// t = VSHR(t, 63);			  // t is 0 (zero) or 1 (not zero)
	// t = VXOR(t, VSET1(1));		  // t = t ^ 1
	return _mm512_reduce_or_epi64(t);

	// return t;
}
