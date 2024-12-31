#include "pair_ifma.h"

void ep2_norm_projc(ep2_t r, const ep2_t p) {
	if (p->coord != BASIC) {
		fp2_t t0, t1;

		fp2_null(t0);
		fp2_null(t1);
		fp2_inv(t1, p->z);

        fp2_mul(r->x, p->x, t1);

        fp2_mul(r->y, p->y, t1);
        fp2_set_dig(r->z, 1);


        fp2_free(t0);
        fp2_free(t1);

	}

	r->coord = BASIC;
}

/* *(1+i) */
void fp2_mul_1i(fp2_t r, fp2_t  a) {
    fp2_t c;
    fp_sub(c[0], a[0], a[1]);
    fp_add(c[1], a[0], a[1]);

    fp_copy(r[0], c[0]);
    fp_copy(r[1], c[1]);
}

void ep_add_projective_mix(ep2_t r, const ep2_t p, const ep2_t q)
{
    fp2_t t0, t1, t2, t3, t4, t5, b3;
    fp2_null(t0);
	fp2_null(t1);
	fp2_null(t2);
	fp2_null(t3);
	fp2_null(t4);
	fp2_null(t5);

	RLC_TRY {
		fp2_new(t0);
		fp2_new(t1);
		fp2_new(t2);
		fp2_new(t3);
		fp2_new(t4);
		fp2_new(t5);

    	fp2_mul(t0, p->x, q->x);
        fp2_mul(t1, p->y, q->y);
        fp2_add(t3, q->x, q->y);
        fp2_add(t4, p->x, p->y);
        fp2_mul(t3, t3, t4);
        fp2_add(t4, t0, t1);
        fp2_sub(t3, t3, t4);

        if (p->coord == BASIC) {
            /* Save 1M + 1m_3b if z1 = 1. */
            fp2_dbl(b3, ep2_curve_get_b());
            fp2_add(b3, b3, ep2_curve_get_b());

            fp2_add(t4, q->y, p->y);
            fp2_add(r->y, q->x, p->x);
            fp2_add(r->z, t1, b3);
            fp2_sub(t1, t1, b3);
		} else {
            /* Cost of 11M + 2m_3b + 13a. */
            fp2_mul(t4, q->y, p->z);
            fp2_add(t4, t4, p->y);
            fp2_mul(r->y, q->x, p->z);
            fp2_add(r->y, r->y, p->x);
            // ep_curve_mul_b3(t2, p->z);
            fp2_mul_1i(t2, p->z);
            fp2_dbl(t5, t2);
            fp2_add(t2, t2, t5);
            fp2_dbl(t2, t2);
            fp2_dbl(t2, t2);
            fp2_add(r->z, t1, t2);
            fp2_sub(t1, t1, t2);
        }

        fp2_dbl(r->x, t0);
        fp2_add(t0, t0, r->x);
        // ep_curve_mul_b3(r->y, r->y);
        fp2_mul_1i(r->y, r->y);
        fp2_dbl(t5, r->y);
        fp2_add(r->y, r->y, t5);
        fp2_dbl(r->y, r->y);
        fp2_dbl(r->y, r->y);
        fp2_mul(r->x, t4, r->y);
        fp2_mul(t2, t3, t1);
        fp2_sub(r->x, t2, r->x);
        fp2_mul(r->y, t0, r->y);
        fp2_mul(t1, t1, r->z);
        fp2_add(r->y, t1, r->y);
        fp2_mul(t0, t0, t3);
        fp2_mul(r->z, r->z, t4);
        fp2_add(r->z, r->z, t0);

        r->coord = PROJC;
    }
    RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp_free(t0);
		fp_free(t1);
		fp_free(t2);
		fp_free(t3);
		fp_free(t4);
		fp_free(t5);
	}
}

void ep2_add_projective_imp(ep2_t r, const ep2_t p, const ep2_t q)
{
//     #if defined(EP_MIXED) || !defined(STRIP)
// 	/* Test if z2 = 1 only if mixed coordinates are turned on. */
	if (q->coord == BASIC) {
		ep_add_projective_mix(r, p, q);
		return;
	}
// #endif
    fp2_t t0, t1, t2, t3, t4, t5;
    fp2_null(t0);
	fp2_null(t1);
	fp2_null(t2);
	fp2_null(t3);
	fp2_null(t4);
	fp2_null(t5);
	fp2_null(t6);

	RLC_TRY {
		fp2_new(t0);
		fp2_new(t1);
		fp2_new(t2);
		fp2_new(t3);
		fp2_new(t4);
		fp2_new(t5);
		fp2_new(t6);

		fp2_mul(t0, p->x, q->x);
		fp2_mul(t1, p->y, q->y);
		fp2_mul(t2, p->z, q->z);
		fp2_add(t3, p->x, p->y);
		fp2_add(t4, q->x, q->y);
		fp2_mul(t3, t3, t4);
		fp2_add(t4, t0, t1);
		fp2_sub(t3, t3, t4);

        fp2_add(t4, p->y, p->z);
        fp2_add(t5, q->y, q->z);
        fp2_mul(t4, t4, t5);
        fp2_add(t5, t1, t2);
        fp2_sub(t4, t4, t5);
        fp2_add(r->y, q->x, q->z);
        fp2_add(r->x, p->x, p->z);
        fp2_mul(r->x, r->x, r->y);
        fp2_add(r->y, t0, t2);
        fp2_sub(r->y, r->x, r->y);
        fp2_dbl(r->x, t0);
        fp2_add(t0, t0, r->x);
        // ep_curve_mul_b3(t2, t2);
        fp2_mul_1i(t2, t2);
        fp2_dbl(t5, t2);
        fp2_add(t2, t2, t5);
        fp2_dbl(t2, t2);
        fp2_dbl(t2, t2);

        fp2_add(r->z, t1, t2);
        fp2_sub(t1, t1, t2);
        // ep_curve_mul_b3(r->y, r->y);
        fp2_mul_1i(r->y, r->y);
        fp2_dbl(t5, r->y);
        fp2_add(r->y, r->y, t5);
        fp2_dbl(r->y, r->y);
        fp2_dbl(r->y, r->y);

        fp2_mul(r->x, t4, r->y);
        fp2_mul(t2, t3, t1);
        fp2_sub(r->x, t2, r->x);
        fp2_mul(r->y, t0, r->y);
        fp2_mul(t1, t1, r->z);
        fp2_add(r->y, t1, r->y);
        fp2_mul(t0, t0, t3);
        fp2_mul(r->z, r->z, t4);
        fp2_add(r->z, r->z, t0);
		r->coord = PROJC;
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(t0);
		fp2_free(t1);
		fp2_free(t2);
		fp2_free(t3);
		fp2_free(t4);
		fp2_free(t5);
		fp2_free(t6);
    }
}

void ep2_add_projective(ep2_t r, const ep2_t p, const ep2_t q) 
{
    if (ep2_is_infty(p)) {
		ep2_copy(r, q);
		return;
	}

	if (ep2_is_infty(q)) {
		ep2_copy(r, p);
		return;
	}

    ep2_add_projective_imp(r, p, q);

}

void ep2_dbl_projective_imp(ep2_t r, const ep2_t p)
{
    fp2_t t0, t1, t2, t3, t4, t5;

	fp2_null(t0);
	fp2_null(t1);
	fp2_null(t2);
	fp2_null(t3);
	fp2_null(t4);
	fp2_null(t5);

    RLC_TRY {
		fp2_new(t0);
		fp2_new(t1);
		fp2_new(t2);
		fp2_new(t3);
		fp2_new(t4);
		fp2_new(t5);

        fp2_sqr(t0, p->y);
        fp2_mul(t3, p->x, p->y);
        if (p->coord == BASIC) {
            /* Save 1M + 1S + 1m_b3 if z1 = 1. */
            fp2_copy(t1, p->y);
            // fp2_copy(t2, ep2_curve_get_b3());
            fp2_copy(t2, ep2_curve_get_b());
            fp2_dbl(t5, t2);
            fp2_add(t2, t2, t5);
        } else {
            fp2_mul(t1, p->y, p->z);
            fp2_sqr(t2, p->z);
            // ep_curve_mul_b3(t2, t2);
            fp2_mul_1i(t2, t2);
            fp2_dbl(t5, t2);
            fp2_add(t2, t2, t5);
            fp2_dbl(t2, t2);
            fp2_dbl(t2, t2);
        }
        fp2_dbl(r->z, t0);
        fp2_dbl(r->z, r->z);
        fp2_dbl(r->z, r->z);
        fp2_mul(r->x, t2, r->z);
        fp2_add(r->y, t0, t2);
        fp2_mul(r->z, t1, r->z);
        fp2_dbl(t1, t2);
        fp2_add(t2, t1, t2);
        fp2_sub(t0, t0, t2);
        fp2_mul(r->y, t0, r->y);
        fp2_add(r->y, r->x, r->y);
        fp2_mul(r->x, t0, t3);
        fp2_dbl(r->x, r->x);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp_free(t0);
		fp_free(t1);
		fp_free(t2);
		fp_free(t3);
		fp_free(t4);
		fp_free(t5);
	}
}


void ep2_dbl_projective(ep2_t r, const ep2_t p) {
	if (ep2_is_infty(p)) {
		ep2_set_infty(r);
		return;
	}

	ep2_dbl_projective_imp(r, p);
}


void simdPADDaff_G1(uint64_t u1x[][HT_NWORDS], uint64_t u1y[][HT_NWORDS], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], pair52 p3[8])
{
    htpoint a, b, c;
    int i;

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(u1x[7][i], u1x[6][i], u1x[5][i], u1x[4][i], u1x[3][i], u1x[2][i], u1x[1][i], u1x[0][i]);
        a.y[i] = set_vector(u1y[7][i], u1y[6][i], u1y[5][i], u1y[4][i], u1y[3][i], u1y[2][i], u1y[1][i], u1y[0][i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
    }

    gfp_mont_relic2avx_8x1w(a.x, a.x);
    gfp_mont_relic2avx_8x1w(a.y, a.y);
    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.z[i] = VSET1(ht_montR[i]);
        b.z[i] = VSET1(ht_montR[i]);
    }
    PADD_projc_8x1w_z1eqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->g1->x, c.x, 0);
    get_channel_8x1w(p3[1]->g1->x, c.x, 1);
    get_channel_8x1w(p3[2]->g1->x, c.x, 2);
    get_channel_8x1w(p3[3]->g1->x, c.x, 3);
    get_channel_8x1w(p3[4]->g1->x, c.x, 4);
    get_channel_8x1w(p3[5]->g1->x, c.x, 5);
    get_channel_8x1w(p3[6]->g1->x, c.x, 6);
    get_channel_8x1w(p3[7]->g1->x, c.x, 7);

    get_channel_8x1w(p3[0]->g1->y, c.y, 0);
    get_channel_8x1w(p3[1]->g1->y, c.y, 1);
    get_channel_8x1w(p3[2]->g1->y, c.y, 2);
    get_channel_8x1w(p3[3]->g1->y, c.y, 3);
    get_channel_8x1w(p3[4]->g1->y, c.y, 4);
    get_channel_8x1w(p3[5]->g1->y, c.y, 5);
    get_channel_8x1w(p3[6]->g1->y, c.y, 6);
    get_channel_8x1w(p3[7]->g1->y, c.y, 7);

    get_channel_8x1w(p3[0]->g1->z, c.z, 0);
    get_channel_8x1w(p3[1]->g1->z, c.z, 1);
    get_channel_8x1w(p3[2]->g1->z, c.z, 2);
    get_channel_8x1w(p3[3]->g1->z, c.z, 3);
    get_channel_8x1w(p3[4]->g1->z, c.z, 4);
    get_channel_8x1w(p3[5]->g1->z, c.z, 5);
    get_channel_8x1w(p3[6]->g1->z, c.z, 6);
    get_channel_8x1w(p3[7]->g1->z, c.z, 7);
    
}

void simdPADDaff_G2(uint64_t u1x[][2][HT_NWORDS], uint64_t u1y[][2][HT_NWORDS], uint64_t u2x[][2][HT_NWORDS], uint64_t u2y[][2][HT_NWORDS], pair52 p3[8])
{
    ht2point a, b, c;
    int i;

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[0][i] = set_vector(u1x[7][0][i], u1x[6][0][i], u1x[5][0][i], u1x[4][0][i], u1x[3][0][i], u1x[2][0][i], u1x[1][0][i], u1x[0][0][i]);
        a.y[0][i] = set_vector(u1y[7][0][i], u1y[6][0][i], u1y[5][0][i], u1y[4][0][i], u1y[3][0][i], u1y[2][0][i], u1y[1][0][i], u1y[0][0][i]);
        b.x[0][i] = set_vector(u2x[7][0][i], u2x[6][0][i], u2x[5][0][i], u2x[4][0][i], u2x[3][0][i], u2x[2][0][i], u2x[1][0][i], u2x[0][0][i]);
        b.y[0][i] = set_vector(u2y[7][0][i], u2y[6][0][i], u2y[5][0][i], u2y[4][0][i], u2y[3][0][i], u2y[2][0][i], u2y[1][0][i], u2y[0][0][i]);
        a.x[1][i] = set_vector(u1x[7][1][i], u1x[6][1][i], u1x[5][1][i], u1x[4][1][i], u1x[3][1][i], u1x[2][1][i], u1x[1][1][i], u1x[0][1][i]);
        a.y[1][i] = set_vector(u1y[7][1][i], u1y[6][1][i], u1y[5][1][i], u1y[4][1][i], u1y[3][1][i], u1y[2][1][i], u1y[1][1][i], u1y[0][1][i]);
        b.x[1][i] = set_vector(u2x[7][1][i], u2x[6][1][i], u2x[5][1][i], u2x[4][1][i], u2x[3][1][i], u2x[2][1][i], u2x[1][1][i], u2x[0][1][i]);
        b.y[1][i] = set_vector(u2y[7][1][i], u2y[6][1][i], u2y[5][1][i], u2y[4][1][i], u2y[3][1][i], u2y[2][1][i], u2y[1][1][i], u2y[0][1][i]);
    }

    gfp2_mont_relic2avx_8x1w(a.x, a.x);
    gfp2_mont_relic2avx_8x1w(a.y, a.y);
    gfp2_mont_relic2avx_8x1w(b.x, b.x);
    gfp2_mont_relic2avx_8x1w(b.y, b.y);

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.z[0][i] = VSET1(ht_montR[i]);
        b.z[0][i] = VSET1(ht_montR[i]);
        a.z[1][i] = VZERO;
        b.z[1][i] = VZERO;
    }

    PADD_G2_projc_8x1w_z1eqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->g2->x[0], c.x[0], 0);
    get_channel_8x1w(p3[1]->g2->x[0], c.x[0], 1);
    get_channel_8x1w(p3[2]->g2->x[0], c.x[0], 2);
    get_channel_8x1w(p3[3]->g2->x[0], c.x[0], 3);
    get_channel_8x1w(p3[4]->g2->x[0], c.x[0], 4);
    get_channel_8x1w(p3[5]->g2->x[0], c.x[0], 5);
    get_channel_8x1w(p3[6]->g2->x[0], c.x[0], 6);
    get_channel_8x1w(p3[7]->g2->x[0], c.x[0], 7);

    get_channel_8x1w(p3[0]->g2->x[1], c.x[1], 0);
    get_channel_8x1w(p3[1]->g2->x[1], c.x[1], 1);
    get_channel_8x1w(p3[2]->g2->x[1], c.x[1], 2);
    get_channel_8x1w(p3[3]->g2->x[1], c.x[1], 3);
    get_channel_8x1w(p3[4]->g2->x[1], c.x[1], 4);
    get_channel_8x1w(p3[5]->g2->x[1], c.x[1], 5);
    get_channel_8x1w(p3[6]->g2->x[1], c.x[1], 6);
    get_channel_8x1w(p3[7]->g2->x[1], c.x[1], 7);

    get_channel_8x1w(p3[0]->g2->y[0], c.y[0], 0);
    get_channel_8x1w(p3[1]->g2->y[0], c.y[0], 1);
    get_channel_8x1w(p3[2]->g2->y[0], c.y[0], 2);
    get_channel_8x1w(p3[3]->g2->y[0], c.y[0], 3);
    get_channel_8x1w(p3[4]->g2->y[0], c.y[0], 4);
    get_channel_8x1w(p3[5]->g2->y[0], c.y[0], 5);
    get_channel_8x1w(p3[6]->g2->y[0], c.y[0], 6);
    get_channel_8x1w(p3[7]->g2->y[0], c.y[0], 7);

    get_channel_8x1w(p3[0]->g2->y[1], c.y[1], 0);
    get_channel_8x1w(p3[1]->g2->y[1], c.y[1], 1);
    get_channel_8x1w(p3[2]->g2->y[1], c.y[1], 2);
    get_channel_8x1w(p3[3]->g2->y[1], c.y[1], 3);
    get_channel_8x1w(p3[4]->g2->y[1], c.y[1], 4);
    get_channel_8x1w(p3[5]->g2->y[1], c.y[1], 5);
    get_channel_8x1w(p3[6]->g2->y[1], c.y[1], 6);
    get_channel_8x1w(p3[7]->g2->y[1], c.y[1], 7);


    get_channel_8x1w(p3[0]->g2->z[0], c.z[0], 0);
    get_channel_8x1w(p3[1]->g2->z[0], c.z[0], 1);
    get_channel_8x1w(p3[2]->g2->z[0], c.z[0], 2);
    get_channel_8x1w(p3[3]->g2->z[0], c.z[0], 3);
    get_channel_8x1w(p3[4]->g2->z[0], c.z[0], 4);
    get_channel_8x1w(p3[5]->g2->z[0], c.z[0], 5);
    get_channel_8x1w(p3[6]->g2->z[0], c.z[0], 6);
    get_channel_8x1w(p3[7]->g2->z[0], c.z[0], 7);

    get_channel_8x1w(p3[0]->g2->z[1], c.z[1], 0);
    get_channel_8x1w(p3[1]->g2->z[1], c.z[1], 1);
    get_channel_8x1w(p3[2]->g2->z[1], c.z[1], 2);
    get_channel_8x1w(p3[3]->g2->z[1], c.z[1], 3);
    get_channel_8x1w(p3[4]->g2->z[1], c.z[1], 4);
    get_channel_8x1w(p3[5]->g2->z[1], c.z[1], 5);
    get_channel_8x1w(p3[6]->g2->z[1], c.z[1], 6);
    get_channel_8x1w(p3[7]->g2->z[1], c.z[1], 7);
    
}

void simdPADDprj_G1(pair_buf2_st* buf, pair52 p3[8])
{
    htpoint a, b, c;
    int i;

    for(i = 1; i < HT_NWORDS; i++)
    {
        if(buf[i].op1 == NULL)
        {
            buf[i].op1 = buf[0].op1;
            buf[i].op2 = buf[0].op2;
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector((*buf[7].op1)->g1->x[i], (*buf[6].op1)->g1->x[i], (*buf[5].op1)->g1->x[i], (*buf[4].op1)->g1->x[i], (*buf[3].op1)->g1->x[i], (*buf[2].op1)->g1->x[i], (*buf[1].op1)->g1->x[i], (*buf[0].op1)->g1->x[i]);
        a.y[i] = set_vector((*buf[7].op1)->g1->y[i], (*buf[6].op1)->g1->y[i], (*buf[5].op1)->g1->y[i], (*buf[4].op1)->g1->y[i], (*buf[3].op1)->g1->y[i], (*buf[2].op1)->g1->y[i], (*buf[1].op1)->g1->y[i], (*buf[0].op1)->g1->y[i]);
        a.z[i] = set_vector((*buf[7].op1)->g1->z[i], (*buf[6].op1)->g1->z[i], (*buf[5].op1)->g1->z[i], (*buf[4].op1)->g1->z[i], (*buf[3].op1)->g1->z[i], (*buf[2].op1)->g1->z[i], (*buf[1].op1)->g1->z[i], (*buf[0].op1)->g1->z[i]);
        b.x[i] = set_vector((*buf[7].op2)->g1->x[i], (*buf[6].op2)->g1->x[i], (*buf[5].op2)->g1->x[i], (*buf[4].op2)->g1->x[i], (*buf[3].op2)->g1->x[i], (*buf[2].op2)->g1->x[i], (*buf[1].op2)->g1->x[i], (*buf[0].op2)->g1->x[i]);
        b.y[i] = set_vector((*buf[7].op2)->g1->y[i], (*buf[6].op2)->g1->y[i], (*buf[5].op2)->g1->y[i], (*buf[4].op2)->g1->y[i], (*buf[3].op2)->g1->y[i], (*buf[2].op2)->g1->y[i], (*buf[1].op2)->g1->y[i], (*buf[0].op2)->g1->y[i]);
        b.z[i] = set_vector((*buf[7].op2)->g1->z[i], (*buf[6].op2)->g1->z[i], (*buf[5].op2)->g1->z[i], (*buf[4].op2)->g1->z[i], (*buf[3].op2)->g1->z[i], (*buf[2].op2)->g1->z[i], (*buf[1].op2)->g1->z[i], (*buf[0].op2)->g1->z[i]);
    }


    PADD_projc_8x1w_z1neqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->g1->x, c.x, 0);
    get_channel_8x1w(p3[1]->g1->x, c.x, 1);
    get_channel_8x1w(p3[2]->g1->x, c.x, 2);
    get_channel_8x1w(p3[3]->g1->x, c.x, 3);
    get_channel_8x1w(p3[4]->g1->x, c.x, 4);
    get_channel_8x1w(p3[5]->g1->x, c.x, 5);
    get_channel_8x1w(p3[6]->g1->x, c.x, 6);
    get_channel_8x1w(p3[7]->g1->x, c.x, 7);

    get_channel_8x1w(p3[0]->g1->y, c.y, 0);
    get_channel_8x1w(p3[1]->g1->y, c.y, 1);
    get_channel_8x1w(p3[2]->g1->y, c.y, 2);
    get_channel_8x1w(p3[3]->g1->y, c.y, 3);
    get_channel_8x1w(p3[4]->g1->y, c.y, 4);
    get_channel_8x1w(p3[5]->g1->y, c.y, 5);
    get_channel_8x1w(p3[6]->g1->y, c.y, 6);
    get_channel_8x1w(p3[7]->g1->y, c.y, 7);

    get_channel_8x1w(p3[0]->g1->z, c.z, 0);
    get_channel_8x1w(p3[1]->g1->z, c.z, 1);
    get_channel_8x1w(p3[2]->g1->z, c.z, 2);
    get_channel_8x1w(p3[3]->g1->z, c.z, 3);
    get_channel_8x1w(p3[4]->g1->z, c.z, 4);
    get_channel_8x1w(p3[5]->g1->z, c.z, 5);
    get_channel_8x1w(p3[6]->g1->z, c.z, 6);
    get_channel_8x1w(p3[7]->g1->z, c.z, 7);
    
}

void simdPADDprj_G2(pair_buf2_st* buf, pair52 p3[8])
{
    ht2point a, b, c;
    int i;


    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[0][i] = set_vector((*buf[7].op1)->g2->x[0][i], (*buf[6].op1)->g2->x[0][i], (*buf[5].op1)->g2->x[0][i], (*buf[4].op1)->g2->x[0][i], (*buf[3].op1)->g2->x[0][i], (*buf[2].op1)->g2->x[0][i], (*buf[1].op1)->g2->x[0][i], (*buf[0].op1)->g2->x[0][i]);
        a.x[1][i] = set_vector((*buf[7].op1)->g2->x[1][i], (*buf[6].op1)->g2->x[1][i], (*buf[5].op1)->g2->x[1][i], (*buf[4].op1)->g2->x[1][i], (*buf[3].op1)->g2->x[1][i], (*buf[2].op1)->g2->x[1][i], (*buf[1].op1)->g2->x[1][i], (*buf[0].op1)->g2->x[1][i]);
        a.y[0][i] = set_vector((*buf[7].op1)->g2->y[0][i], (*buf[6].op1)->g2->y[0][i], (*buf[5].op1)->g2->y[0][i], (*buf[4].op1)->g2->y[0][i], (*buf[3].op1)->g2->y[0][i], (*buf[2].op1)->g2->y[0][i], (*buf[1].op1)->g2->y[0][i], (*buf[0].op1)->g2->y[0][i]);
        a.y[1][i] = set_vector((*buf[7].op1)->g2->y[1][i], (*buf[6].op1)->g2->y[1][i], (*buf[5].op1)->g2->y[1][i], (*buf[4].op1)->g2->y[1][i], (*buf[3].op1)->g2->y[1][i], (*buf[2].op1)->g2->y[1][i], (*buf[1].op1)->g2->y[1][i], (*buf[0].op1)->g2->y[1][i]);
        a.z[0][i] = set_vector((*buf[7].op1)->g2->z[0][i], (*buf[6].op1)->g2->z[0][i], (*buf[5].op1)->g2->z[0][i], (*buf[4].op1)->g2->z[0][i], (*buf[3].op1)->g2->z[0][i], (*buf[2].op1)->g2->z[0][i], (*buf[1].op1)->g2->z[0][i], (*buf[0].op1)->g2->z[0][i]);
        a.z[1][i] = set_vector((*buf[7].op1)->g2->z[1][i], (*buf[6].op1)->g2->z[1][i], (*buf[5].op1)->g2->z[1][i], (*buf[4].op1)->g2->z[1][i], (*buf[3].op1)->g2->z[1][i], (*buf[2].op1)->g2->z[1][i], (*buf[1].op1)->g2->z[1][i], (*buf[0].op1)->g2->z[1][i]);
        b.x[0][i] = set_vector((*buf[7].op2)->g2->x[0][i], (*buf[6].op2)->g2->x[0][i], (*buf[5].op2)->g2->x[0][i], (*buf[4].op2)->g2->x[0][i], (*buf[3].op2)->g2->x[0][i], (*buf[2].op2)->g2->x[0][i], (*buf[1].op2)->g2->x[0][i], (*buf[0].op2)->g2->x[0][i]);
        b.x[1][i] = set_vector((*buf[7].op2)->g2->x[1][i], (*buf[6].op2)->g2->x[1][i], (*buf[5].op2)->g2->x[1][i], (*buf[4].op2)->g2->x[1][i], (*buf[3].op2)->g2->x[1][i], (*buf[2].op2)->g2->x[1][i], (*buf[1].op2)->g2->x[1][i], (*buf[0].op2)->g2->x[1][i]);
        b.y[0][i] = set_vector((*buf[7].op2)->g2->y[0][i], (*buf[6].op2)->g2->y[0][i], (*buf[5].op2)->g2->y[0][i], (*buf[4].op2)->g2->y[0][i], (*buf[3].op2)->g2->y[0][i], (*buf[2].op2)->g2->y[0][i], (*buf[1].op2)->g2->y[0][i], (*buf[0].op2)->g2->y[0][i]);
        b.y[1][i] = set_vector((*buf[7].op2)->g2->y[1][i], (*buf[6].op2)->g2->y[1][i], (*buf[5].op2)->g2->y[1][i], (*buf[4].op2)->g2->y[1][i], (*buf[3].op2)->g2->y[1][i], (*buf[2].op2)->g2->y[1][i], (*buf[1].op2)->g2->y[1][i], (*buf[0].op2)->g2->y[1][i]);
        b.z[0][i] = set_vector((*buf[7].op2)->g2->z[0][i], (*buf[6].op2)->g2->z[0][i], (*buf[5].op2)->g2->z[0][i], (*buf[4].op2)->g2->z[0][i], (*buf[3].op2)->g2->z[0][i], (*buf[2].op2)->g2->z[0][i], (*buf[1].op2)->g2->z[0][i], (*buf[0].op2)->g2->z[0][i]);
        b.z[1][i] = set_vector((*buf[7].op2)->g2->z[1][i], (*buf[6].op2)->g2->z[1][i], (*buf[5].op2)->g2->z[1][i], (*buf[4].op2)->g2->z[1][i], (*buf[3].op2)->g2->z[1][i], (*buf[2].op2)->g2->z[1][i], (*buf[1].op2)->g2->z[1][i], (*buf[0].op2)->g2->z[1][i]);
    }


    PADD_G2_projc_8x1w_z1neqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->g2->x[0], c.x[0], 0);
    get_channel_8x1w(p3[1]->g2->x[0], c.x[0], 1);
    get_channel_8x1w(p3[2]->g2->x[0], c.x[0], 2);
    get_channel_8x1w(p3[3]->g2->x[0], c.x[0], 3);
    get_channel_8x1w(p3[4]->g2->x[0], c.x[0], 4);
    get_channel_8x1w(p3[5]->g2->x[0], c.x[0], 5);
    get_channel_8x1w(p3[6]->g2->x[0], c.x[0], 6);
    get_channel_8x1w(p3[7]->g2->x[0], c.x[0], 7);

    get_channel_8x1w(p3[0]->g2->x[1], c.x[1], 0);
    get_channel_8x1w(p3[1]->g2->x[1], c.x[1], 1);
    get_channel_8x1w(p3[2]->g2->x[1], c.x[1], 2);
    get_channel_8x1w(p3[3]->g2->x[1], c.x[1], 3);
    get_channel_8x1w(p3[4]->g2->x[1], c.x[1], 4);
    get_channel_8x1w(p3[5]->g2->x[1], c.x[1], 5);
    get_channel_8x1w(p3[6]->g2->x[1], c.x[1], 6);
    get_channel_8x1w(p3[7]->g2->x[1], c.x[1], 7);

    get_channel_8x1w(p3[0]->g2->y[0], c.y[0], 0);
    get_channel_8x1w(p3[1]->g2->y[0], c.y[0], 1);
    get_channel_8x1w(p3[2]->g2->y[0], c.y[0], 2);
    get_channel_8x1w(p3[3]->g2->y[0], c.y[0], 3);
    get_channel_8x1w(p3[4]->g2->y[0], c.y[0], 4);
    get_channel_8x1w(p3[5]->g2->y[0], c.y[0], 5);
    get_channel_8x1w(p3[6]->g2->y[0], c.y[0], 6);
    get_channel_8x1w(p3[7]->g2->y[0], c.y[0], 7);

    get_channel_8x1w(p3[0]->g2->y[1], c.y[1], 0);
    get_channel_8x1w(p3[1]->g2->y[1], c.y[1], 1);
    get_channel_8x1w(p3[2]->g2->y[1], c.y[1], 2);
    get_channel_8x1w(p3[3]->g2->y[1], c.y[1], 3);
    get_channel_8x1w(p3[4]->g2->y[1], c.y[1], 4);
    get_channel_8x1w(p3[5]->g2->y[1], c.y[1], 5);
    get_channel_8x1w(p3[6]->g2->y[1], c.y[1], 6);
    get_channel_8x1w(p3[7]->g2->y[1], c.y[1], 7);


    get_channel_8x1w(p3[0]->g2->z[0], c.z[0], 0);
    get_channel_8x1w(p3[1]->g2->z[0], c.z[0], 1);
    get_channel_8x1w(p3[2]->g2->z[0], c.z[0], 2);
    get_channel_8x1w(p3[3]->g2->z[0], c.z[0], 3);
    get_channel_8x1w(p3[4]->g2->z[0], c.z[0], 4);
    get_channel_8x1w(p3[5]->g2->z[0], c.z[0], 5);
    get_channel_8x1w(p3[6]->g2->z[0], c.z[0], 6);
    get_channel_8x1w(p3[7]->g2->z[0], c.z[0], 7);

    get_channel_8x1w(p3[0]->g2->z[1], c.z[1], 0);
    get_channel_8x1w(p3[1]->g2->z[1], c.z[1], 1);
    get_channel_8x1w(p3[2]->g2->z[1], c.z[1], 2);
    get_channel_8x1w(p3[3]->g2->z[1], c.z[1], 3);
    get_channel_8x1w(p3[4]->g2->z[1], c.z[1], 4);
    get_channel_8x1w(p3[5]->g2->z[1], c.z[1], 5);
    get_channel_8x1w(p3[6]->g2->z[1], c.z[1], 6);
    get_channel_8x1w(p3[7]->g2->z[1], c.z[1], 7);
    
}

void simdPADDmix_G1(pair52_st *bucket[AVX_WAY], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], int num)
{
    htpoint a, b, c;
    int i;

    for(i = num; i < HT_NWORDS; i++)
    {
        if(bucket[i] == NULL)
        {
            bucket[i] = bucket[0];
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(bucket[7]->g1->x[i], bucket[6]->g1->x[i], bucket[5]->g1->x[i], bucket[4]->g1->x[i], bucket[3]->g1->x[i], bucket[2]->g1->x[i], bucket[1]->g1->x[i], bucket[0]->g1->x[i]);
        a.y[i] = set_vector(bucket[7]->g1->y[i], bucket[6]->g1->y[i], bucket[5]->g1->y[i], bucket[4]->g1->y[i], bucket[3]->g1->y[i], bucket[2]->g1->y[i], bucket[1]->g1->y[i], bucket[0]->g1->y[i]);
        a.z[i] = set_vector(bucket[7]->g1->z[i], bucket[6]->g1->z[i], bucket[5]->g1->z[i], bucket[4]->g1->z[i], bucket[3]->g1->z[i], bucket[2]->g1->z[i], bucket[1]->g1->z[i], bucket[0]->g1->z[i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
        b.z[i] = VSET1(ht_montR[i]);
   }

    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);
    PADD_projc_8x1w_mix(&c, &a, &b);

    for (i = 0; i < num ; i++)
    {
        get_channel_8x1w(bucket[i]->g1->x, c.x, i);
        get_channel_8x1w(bucket[i]->g1->y, c.y, i);
        get_channel_8x1w(bucket[i]->g1->z, c.z, i);
    }
    
}

void simdPADDmix_G2(pair52_st *bucket[AVX_WAY], uint64_t u2x[][2][HT_NWORDS], uint64_t u2y[][2][HT_NWORDS], int num)
{
    ht2point a, b, c;
    int i;


    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[0][i] = set_vector(bucket[7]->g2->x[0][i], bucket[6]->g2->x[0][i], bucket[5]->g2->x[0][i], bucket[4]->g2->x[0][i], bucket[3]->g2->x[0][i], bucket[2]->g2->x[0][i], bucket[1]->g2->x[0][i], bucket[0]->g2->x[0][i]);
        a.x[1][i] = set_vector(bucket[7]->g2->x[1][i], bucket[6]->g2->x[1][i], bucket[5]->g2->x[1][i], bucket[4]->g2->x[1][i], bucket[3]->g2->x[1][i], bucket[2]->g2->x[1][i], bucket[1]->g2->x[1][i], bucket[0]->g2->x[1][i]);
        a.y[0][i] = set_vector(bucket[7]->g2->y[0][i], bucket[6]->g2->y[0][i], bucket[5]->g2->y[0][i], bucket[4]->g2->y[0][i], bucket[3]->g2->y[0][i], bucket[2]->g2->y[0][i], bucket[1]->g2->y[0][i], bucket[0]->g2->y[0][i]);
        a.y[1][i] = set_vector(bucket[7]->g2->y[1][i], bucket[6]->g2->y[1][i], bucket[5]->g2->y[1][i], bucket[4]->g2->y[1][i], bucket[3]->g2->y[1][i], bucket[2]->g2->y[1][i], bucket[1]->g2->y[1][i], bucket[0]->g2->y[1][i]);
        a.z[0][i] = set_vector(bucket[7]->g2->z[0][i], bucket[6]->g2->z[0][i], bucket[5]->g2->z[0][i], bucket[4]->g2->z[0][i], bucket[3]->g2->z[0][i], bucket[2]->g2->z[0][i], bucket[1]->g2->z[0][i], bucket[0]->g2->z[0][i]);
        a.z[1][i] = set_vector(bucket[7]->g2->z[1][i], bucket[6]->g2->z[1][i], bucket[5]->g2->z[1][i], bucket[4]->g2->z[1][i], bucket[3]->g2->z[1][i], bucket[2]->g2->z[1][i], bucket[1]->g2->z[1][i], bucket[0]->g2->z[1][i]);
        b.x[0][i] = set_vector(u2x[7][0][i], u2x[6][0][i], u2x[5][0][i], u2x[4][0][i], u2x[3][0][i], u2x[2][0][i], u2x[1][0][i], u2x[0][0][i]);
        b.y[0][i] = set_vector(u2y[7][0][i], u2y[6][0][i], u2y[5][0][i], u2y[4][0][i], u2y[3][0][i], u2y[2][0][i], u2y[1][0][i], u2y[0][0][i]);
        b.x[1][i] = set_vector(u2x[7][1][i], u2x[6][1][i], u2x[5][1][i], u2x[4][1][i], u2x[3][1][i], u2x[2][1][i], u2x[1][1][i], u2x[0][1][i]);
        b.y[1][i] = set_vector(u2y[7][1][i], u2y[6][1][i], u2y[5][1][i], u2y[4][1][i], u2y[3][1][i], u2y[2][1][i], u2y[1][1][i], u2y[0][1][i]);
        b.z[0][i] = VSET1(ht_montR[i]);
        b.z[1][i] = VZERO;
   }

    gfp2_mont_relic2avx_8x1w(b.x, b.x);
    gfp2_mont_relic2avx_8x1w(b.y, b.y);
    PADD_G2_projc_8x1w_mix(&c, &a, &b);

    for (i = 0; i < num ; i++)
    {
        get_channel_8x1w(bucket[i]->g2->x[0], c.x[0], i);
        get_channel_8x1w(bucket[i]->g2->x[1], c.x[1], i);
        get_channel_8x1w(bucket[i]->g2->y[0], c.y[0], i);
        get_channel_8x1w(bucket[i]->g2->y[1], c.y[1], i);
        get_channel_8x1w(bucket[i]->g2->z[0], c.z[0], i);
        get_channel_8x1w(bucket[i]->g2->z[1], c.z[1], i);
    }
    
}



void avx2_pair_withIFMA(pair_msm_ctx *ctx, int eff_times)
{
	int i;
    pair_buf2_st *item;
    pair52 avx2_out[AVX_WAY];

    simdPADDprj_G1(ctx->buf2, avx2_out);
    simdPADDprj_G2(ctx->buf2, avx2_out);

    for (i = 0; i < eff_times; i++)
    {
        item = &(ctx->buf2[i]);
        if(item->flag == 0)
		{
			pair_queue_in(&(ctx->T), avx2_out[i], item->id);
		}
		else
		{
			// ep_copy(ctx->bucket[item->id].B, avx2_out[i]);
            pair52_copy(*item->op2, avx2_out[i]);
		}
    }

    ctx->cnt2 = 0;
}

buf2_st* check_eq_pair_id(pair_msm_ctx *ctx, int id) {
    int i = 0;
    while (i < ctx->cnt2) {
        if ((ctx->buf2[i].id == id) && (ctx->buf2[i].flag == 1)) {
            return &(ctx->buf2[i]);
        }
        i++;
    }
    return NULL;
}

void pair_state_trans_T_withIFMA(pair_msm_ctx *ctx)
{
	pair_buf2_st *same_id;
	int bucket_id, i;

	while (!pair_queue_is_empty(&(ctx->T))) 
	{
		i = ctx->T.front;
        bucket_id = ctx->T.pbuf[i].id;
		same_id = check_eq_pair_id(ctx, bucket_id);
		if(same_id == NULL)
		{	

            /* write into bucket*/
            ctx->buf2[ctx->cnt2].id = bucket_id;
            ctx->buf2[ctx->cnt2].op1 = &(ctx->T.pbuf[i].point);
            ctx->buf2[ctx->cnt2].op2 = &(ctx->bucket[bucket_id].B);
            ctx->buf2[ctx->cnt2].flag = 1;
            ctx->cnt2++;

		}
		else
		{
			/* write into queue*/
	    	same_id->op2 = &(ctx->T.pbuf[i].point);
			same_id->flag = 0;
		}

		if(ctx->cnt2 == AVX_WAY)
		{
			avx2_pair_withIFMA(ctx, AVX_WAY);
		}

		pair_queue_out(&(ctx->T), 1);
    }
}

void avx1_pair_withIFMA(pair_msm_ctx *ctx, int eff_times)
{
	int i;
    pair_buf1_st *dblp;
    pair *p1, *p2;
    pair52 avx1_out[AVX_WAY];
    uint64_t u1x[AVX_WAY][HT_NWORDS], u1y[AVX_WAY][HT_NWORDS], u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    uint64_t v1x[AVX_WAY][2][HT_NWORDS], v1y[AVX_WAY][2][HT_NWORDS], v2x[AVX_WAY][2][HT_NWORDS], v2y[AVX_WAY][2][HT_NWORDS];
 
	/* avx1*/
	for(i = 0; i < eff_times; i++)
	{
		dblp = &ctx->buf1[i];

        p1 = (pair *)(dblp->P0);
        p2 = (pair *)(dblp->P1);

        /* G1*/
        mpi_conv_64to52(u1x[i], (*p1)->g1->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u1y[i], (*p1)->g1->y, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2x[i], (*p2)->g1->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2y[i], (*p2)->g1->y, HT_NWORDS, 6); // convert to radix-52

        /* G2*/
        mpi_conv_64to52(v1x[i][0], (*p1)->g2->x[0], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v1y[i][0], (*p1)->g2->y[0], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v2x[i][0], (*p2)->g2->x[0], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v2y[i][0], (*p2)->g2->y[0], HT_NWORDS, 6); // convert to radix-52

        mpi_conv_64to52(v1x[i][1], (*p1)->g2->x[1], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v1y[i][1], (*p1)->g2->y[1], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v2x[i][1], (*p2)->g2->x[1], HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(v2y[i][1], (*p2)->g2->y[1], HT_NWORDS, 6); // convert to radix-52
	}
    simdPADDaff_G1(u1x, u1y, u2x, u2y, avx1_out);
    simdPADDaff_G2(v1x, v1y, v2x, v2y, avx1_out);

    for(i = 0; i < eff_times; i++)
    {
        dblp = &ctx->buf1[i];
        if(point52_is_infty(ctx->bucket[dblp->id].B->g1))
        {
            pair52_copy(ctx->bucket[dblp->id].B, avx1_out[i]);
        }
        else
        {
            pair_queue_in(&(ctx->T), avx1_out[i], dblp->id);
        }
    }

	ctx->cnt1 = 0;

	/* deal with the avx1_out, change state in T*/
	pair_state_trans_T_withIFMA(ctx);
}


void put_pair_into_bucket_withIFMA(const pair *P, const uint32_t b_id, int p_id, pair_msm_ctx *ctx)
{
	/* judge if the bucket wait is empty*/
	if(ctx->bucket[b_id].wait == NULL)
	{
		// printf("%d ", p_id);
		ctx->bucket[b_id].wait = &(P[p_id]);
	}
	else
	{
		ctx->buf1[ctx->cnt1].id = b_id;
		ctx->buf1[ctx->cnt1].P0 = ctx->bucket[b_id].wait;
		ctx->buf1[ctx->cnt1].P1 = &(P[p_id]);
		ctx->cnt1++;
		ctx->bucket[b_id].wait = NULL;

		/* judge if buf1 is full*/
		if(ctx->cnt1 == AVX_WAY)
		{
			/* avx1_in is full, then do avx1*/
			avx1_pair_withIFMA(ctx, AVX_WAY);
		}
	}
   
}

void pair_mont64_to_mont52(pair52 q, pair p)
{
    ep_t tmp;
    ep2_t tmp2;
    fp_t mont_a;
    mont_a[0] = 0x44f6480ea8e9b9af;
    mont_a[1] = 0xa96f7d65766c8fe4;
    mont_a[2] = 0xe82efd4228b540fe;
    mont_a[3] = 0x6723e5f0ade53b2e;
    mont_a[4] = 0x25ff6eb6fdd4230a;
    mont_a[5] = 0x14c8ee06ef23c24a;

    ep_copy(tmp, p->g1);
    fp_mul(tmp->x, tmp->x, mont_a);
    fp_mul(tmp->y, tmp->y, mont_a);
    fp_mul(tmp->z, tmp->z, mont_a);
    mpi_conv_64to52(q->g1->x, tmp->x, HT_NWORDS, 6);
    mpi_conv_64to52(q->g1->y, tmp->y, HT_NWORDS, 6);
    mpi_conv_64to52(q->g1->z, tmp->z, HT_NWORDS, 6);

    ep2_copy(tmp2, p->g2);
    fp_mul(tmp2->x[0], tmp2->x[0], mont_a);
    fp_mul(tmp2->y[0], tmp2->y[0], mont_a);
    fp_mul(tmp2->z[0], tmp2->z[0], mont_a);
    fp_mul(tmp2->x[1], tmp2->x[1], mont_a);
    fp_mul(tmp2->y[1], tmp2->y[1], mont_a);
    fp_mul(tmp2->z[1], tmp2->z[1], mont_a);
    mpi_conv_64to52(q->g2->x[0], tmp2->x[0], HT_NWORDS, 6);
    mpi_conv_64to52(q->g2->y[0], tmp2->y[0], HT_NWORDS, 6);
    mpi_conv_64to52(q->g2->z[0], tmp2->z[0], HT_NWORDS, 6);
    mpi_conv_64to52(q->g2->x[1], tmp2->x[1], HT_NWORDS, 6);
    mpi_conv_64to52(q->g2->y[1], tmp2->y[1], HT_NWORDS, 6);
    mpi_conv_64to52(q->g2->z[1], tmp2->z[1], HT_NWORDS, 6);
}

void pair_mont52_to_mont64(pair52 q, pair p)
{
    ep_t tmp;
    ep2_t tmp2;
    fp_t mont_b;

    mont_b[0] = 0x0;
    mont_b[1] = 0x0;
    mont_b[2] = 0x0;
    mont_b[3] = 0x0;
    mont_b[4] = 0x0;
    mont_b[5] = 0x100000000;
    ep_set_infty(tmp);
    ep2_set_infty(tmp2);

    mpi_conv_52to64(tmp2->x[0], q->g2->x[0], 6, HT_NWORDS);
    mpi_conv_52to64(tmp2->y[0], q->g2->y[0], 6, HT_NWORDS);
    mpi_conv_52to64(tmp2->z[0], q->g2->z[0], 6, HT_NWORDS);
    mpi_conv_52to64(tmp2->x[1], q->g2->x[1], 6, HT_NWORDS);
    mpi_conv_52to64(tmp2->y[1], q->g2->y[1], 6, HT_NWORDS);
    mpi_conv_52to64(tmp2->z[1], q->g2->z[1], 6, HT_NWORDS);
    fp_mul(p->g2->x[0], tmp2->x[0], mont_b);
    fp_mul(p->g2->x[1], tmp2->x[1], mont_b);
    fp_mul(p->g2->y[0], tmp2->y[0], mont_b);
    fp_mul(p->g2->y[1], tmp2->y[1], mont_b);
    fp_mul(p->g2->z[0], tmp2->z[0], mont_b);
    fp_mul(p->g2->z[1], tmp2->z[1], mont_b);
    
    mpi_conv_52to64(tmp->x, q->g1->x, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->y, q->g1->y, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->z, q->g1->z, 6, HT_NWORDS);
    fp_mul(p->g1->x, tmp->x, mont_b);
    fp_mul(p->g1->y, tmp->y, mont_b);
    fp_mul(p->g1->z, tmp->z, mont_b);
}

void pippenger_pair_ifma(const bn_t *k, const pair *P, int num, pair result, pair *M, pair_msm_ctx *ctx)
{
    int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

	pair Gu[NWINS], Tu[NWINS];
	uint32_t b;
	int i, a;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    uint64_t v2x[AVX_WAY][2][HT_NWORDS], v2y[AVX_WAY][2][HT_NWORDS];
    pair tmp;
    
    for (i = 0; i < NWINS; i++) {
		// printf("第%d列 \n", i);
		init_pair_ctx(ctx);
		for (int j = 0; j < num; j++) {
            b = get_wval(k[j], i, WBITS);
			if (b != 0) {
				// ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
				put_pair_into_bucket_withIFMA(P, (uint32_t)b - 1, j, ctx);
			}
		}

		/* deal with the remains in buf*/
		/* remain buf1*/
		if(ctx->cnt1 != 0)
		{
			avx1_pair_withIFMA(ctx, ctx->cnt1);
		}
		/* remain buf2*/
		while((!pair_queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
		{
			avx2_pair_withIFMA(ctx, ctx->cnt2);
			pair_state_trans_T_withIFMA(ctx);
		}

		/* remain buf0*/
		for(a = 0; a < BUCKETNUM; a++)
		{
            if (ctx->bucket[a].wait != NULL)
            {
                // ep_add_projc(ctx->bucket[ctx->dblpoint[a].id].B, ctx->bucket[ctx->dblpoint[a].id].B, *(ep_t *)(ctx->dblpoint[a].P0));
                if(point52_is_infty(ctx->bucket[a].B->g1))
                {
                    pair_mont64_to_mont52(ctx->bucket[a].B, *(ctx->bucket[a].wait));
                }
                else
                {
                    ctx->avx3_out[ctx->avx3_cnt] = &(ctx->bucket[a].B);
                    mpi_conv_64to52(u2x[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->g1->x, HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(u2y[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->g1->y, HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(v2x[ctx->avx3_cnt][0], (*(ctx->bucket[a].wait))->g2->x[0], HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(v2y[ctx->avx3_cnt][0], (*(ctx->bucket[a].wait))->g2->y[0], HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(v2x[ctx->avx3_cnt][1], (*(ctx->bucket[a].wait))->g2->x[1], HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(v2y[ctx->avx3_cnt][1], (*(ctx->bucket[a].wait))->g2->y[1], HT_NWORDS, 6); // convert to radix-52
                    // mpi_conv_64to52(u2z[ctx->avx3_cnt], (*(ep_t *)(ctx->dblpoint[a].P0))->z, HT_NWORDS, 6); // convert to radix-52
                    ctx->avx3_cnt++;

                    if(ctx->avx3_cnt == AVX_WAY)
                    {
                        simdPADDmix_G1(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
                        simdPADDmix_G2(ctx->avx3_out, v2x, v2y, ctx->avx3_cnt);
                        ctx->avx3_cnt = 0;
                    }
                }
            }
        }
        if(ctx->avx3_cnt)
        {
            simdPADDmix_G1(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
            simdPADDmix_G2(ctx->avx3_out, v2x, v2y, ctx->avx3_cnt);
        }

        pair_mont52_to_mont64(ctx->bucket[BUCKETNUM - 1].B, M[0]);
		// ep_copy(M[0]->g1, ctx->bucket[BUCKETNUM - 1].B->g1);
		// ep2_copy(M[0]->g2, ctx->bucket[BUCKETNUM - 1].B->g2);
        M[0]->g1->coord = PROJC;
        M[0]->g2->coord = PROJC;
        tmp->g1->coord = PROJC;
        tmp->g2->coord = PROJC;
		ep_copy(Gu[i]->g1, M[0]->g1);
		ep2_copy(Gu[i]->g2, M[0]->g2);

		for (int j = 1; j < BUCKETNUM; j++) {
            pair_mont52_to_mont64(ctx->bucket[BUCKETNUM - 1 - j].B, tmp);
			ep_add_projc(M[j]->g1, tmp->g1, M[j - 1]->g1);
			ep2_add_projective(M[j]->g2, tmp->g2, M[j - 1]->g2);
			ep_add_projc(Gu[i]->g1, Gu[i]->g1, M[j]->g1);
			ep2_add_projective(Gu[i]->g2, Gu[i]->g2, M[j]->g2);
		}

        // printf("Gu[%d]:\n", i);
		// ep2_norm_projc(Gu[i]->g2,Gu[i]->g2);
		// ep2_print(Gu[i]->g2);
	}

	ep_copy(Tu[0]->g1, Gu[NWINS - 1]->g1);
    ep2_copy(Tu[0]->g2, Gu[NWINS - 1]->g2);
	for (int i = 1; i < NWINS; i++) {
		for (int j = 0; j < WBITS; j++) {
			ep_dbl_projc(Tu[i - 1]->g1, Tu[i - 1]->g1);
            ep2_dbl_projective(Tu[i - 1]->g2, Tu[i - 1]->g2);
		}
		ep_add_projc(Tu[i]->g1, Gu[NWINS - 1 - i]->g1, Tu[i - 1]->g1);
        ep2_add_projective(Tu[i]->g2, Gu[NWINS - 1 - i]->g2, Tu[i - 1]->g2);
	}
	ep_copy(result->g1, Tu[NWINS - 1]->g1);
    ep2_copy(result->g2, Tu[NWINS - 1]->g2);
    ep2_norm_projc(result->g2, result->g2);
}