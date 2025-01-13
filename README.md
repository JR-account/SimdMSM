# SimdMSM

This source code is an efficient implemetation of MSM and zkSNARK using AVX-512IFMA. It is the artifact of the paper **SimdMSM: SIMD-accelerated Multi-Scalar Multiplication Framework for zkSNARKs** accepted to TCHES 2025. 


## Overview

There are three subfolders included in this repository:
- `AVX-MSM` : the MSM implementation instantiated with AVX512-IFMA engine based on [the RELIC library](https://github.com/relic-toolkit/relic). The specific implementation code can be found in the `AVX-MSM/demo/381/` directory. The AVX512-IFMA engine implementation is based on [Cheng et al.’s work](https://github.com/ulhaocheng/avxcsidh?tab=readme-ov-file).
- `AVX-ZK` : integrating AVX-MSM implementation into [the libsnark library](https://github.com/scipr-lab/libsnark). The part of `r1cs_gg_ppzksnark`, commonly known as the famous Groth16 protocol, is changed to using new AVX-MSM.
- `jsnark` : a tool for evaluating the performance of AVX-ZK under different real-world workloads.

## Requirement

- Ubuntu 22.04.4
- gcc version 11.4.0 
- cmake 3.22.1
- support AVX512-IFMA instruction sets

## Build

### AVX-MSM

 Target the `SimdMSM` library.

```shell
$ cd  AVX-MSM
$ cd demo/381/
$ make lib
```

Run AVX-MSM. The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pip_ifma.c`.

```shell
$ make ifma
$ ./buile/test_pip_ifma
```

Run AVX-pair-MSM. The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pair_ifma.c`.

```shell
$ make pair_ifma
$ ./buile/test_pair_ifma
```

Run AVX-MSM(muti-threads). The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pip_threads.c`.

```shell
$ make ifma
$ ./buile/test_pip_threads
```

Run AVX-pair-MSM(muti-threads). The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pair_threads.c`.

```shell
$ make pair_ifma
$ ./buile/test_pair_threads
```

Generate static link library `libmsm.a`.

```shell
$ make msm

```

### AVX-ZK

Cmake and create the Makefile:

```shell
$ cd AVX-ZK
$ mkdir build && cd build && cmake ..
```

Copy the `libmsm.a` and `librelic_s.a`.to libsnark/build/depends/libff/libff.

Then, to compile the library, run this within the `build` directory:

```shell
$ make
```

Run the profiling of AVX-ZK.

```shell
$ make profile_r1cs_gg_ppzksnark
$ ./libsnark/profile_r1cs_gg_ppzksnark 1000  100 bytes
```

### Running and Testing AVX-ZK by JsnarkCircuitBuilder

To compile the JsnarkCircuitBuilder project via command line, from the jsnark directory:

```shell
$ cd jsnark
$ cd JsnarkCircuitBuilder
$ mkdir -p bin
$ javac -d bin -cp /usr/share/java/junit4.jar:bcprov-jdk15on-159.jar  $(find ./src/* | grep ".java$")
```



Run AES.

```shell
$ java -cp bin examples.generators.blockciphers.AES128CipherCircuitGenerator
```

Run SHA-256.

```shell
$ java -cp bin examples.generators.hash.SHA2CircuitGenerator
```



Run RSAEnc.

```shell
$ java -cp bin examples.generators.rsa.RSAEncryptionCircuitGenerator
```



Run Merkle-Tree.

```shell
$ java -cp bin examples.generators.hash.MerkleTreeMembershipCircuitGenerator
```



Run RSASigVer.

```shell
$ java -cp bin examples.generators.rsa.RSASigVerCircuitGenerator
```

Run Auction.

```shell
$ java -cp bin examples.generators.augmenter.AugmentedAuctionCircuitGenerator
```

