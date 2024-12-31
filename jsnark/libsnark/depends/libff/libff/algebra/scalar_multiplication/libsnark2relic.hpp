#ifndef LIBSNARK2RELIC_HPP_
#define LIBSNARK2RELIC_HPP_

#include <include/relic.h>
#include <libff/algebra/fields/bigint.hpp>
#include <libff/algebra/fields/fp_aux.tcc>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>

namespace libff {

void G1_to_ep(ep_t p, const bls12_381_G1 &q);
void Fr_to_bn(bn_t k, const bls12_381_Fr &vlu);

}
#include <libff/algebra/scalar_multiplication/libsnark2relic.cpp>
#endif // LIBSNARK2RELIC_HPP_