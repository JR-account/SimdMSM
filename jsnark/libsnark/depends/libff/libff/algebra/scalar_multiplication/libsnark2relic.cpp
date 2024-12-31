#ifndef LIBSNARK2RELIC_TCC_
#define LIBSNARK2RELIC_TCC_

#include <libff/algebra/scalar_multiplication/libsnark2relic.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <include/relic.h>


namespace libff {

void G1_to_ep(ep_t p, const bls12_381_G1 &q) {
    // ep_copy(p, *((ep_t*) &q));
    // ep_print(*((ep_t*) &q));
    memcpy(p->x, &(q.X.mont_repr), 384);
    memcpy(p->y, &(q.Y.mont_repr), 384);
    memcpy(p->z, &(q.Z.mont_repr), 384);
    p->coord = BASIC;
}

void Fr_to_bn(bn_t k, const bls12_381_Fr &vlu) {
    // core_init();
    // fp_prime_init();
    // ep_curve_init();
    // fp_param_set(B12_381);
    // ep_param_set(B12_P381);
    // bn_zero(k);
    // memcpy(k->dp, &(vlu.as_bigint()), 4);
    k->dp[0] = vlu.as_bigint().data[0];
    k->dp[1] = vlu.as_bigint().data[1];
    k->dp[2] = vlu.as_bigint().data[2];
    k->dp[3] = vlu.as_bigint().data[3];
    k->used = 4;
    k->sign = RLC_POS;
    // bn_print(k);
}

void ep_to_G1(bls12_381_G1 &q, const ep_t p) {
    ep_copy(*((ep_t*) &q), p);
    q.print();
}

} // libff

#endif