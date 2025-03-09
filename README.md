# SimdMSM

This source code is an efficient implemetation of MSM and zkSNARK using AVX-512IFMA. It is the artifact of the paper [**SimdMSM: SIMD-accelerated Multi-Scalar Multiplication Framework for zkSNARKs**](https://tches.iacr.org/index.php/TCHES/article/view/12061) accepted to [TCHES 2025](https://ches.iacr.org/2025/). 


## Overview

There are three subfolders included in this repository:
- `AVX-MSM` : the MSM implementation instantiated with AVX512-IFMA engine based on [the RELIC library](https://github.com/relic-toolkit/relic). The specific implementation code can be found in the `AVX-MSM/demo/381/` directory. The AVX512-IFMA engine implementation is based on [Cheng et al.’s work](https://github.com/ulhaocheng/avxcsidh?tab=readme-ov-file).
- `AVX-ZK` : integrating AVX-MSM implementation into [the libsnark library](https://github.com/scipr-lab/libsnark). The part of `r1cs_gg_ppzksnark`, commonly known as the famous Groth16 protocol, is changed to using new AVX-MSM.
- `jsnark` : a tool for evaluating the performance of AVX-ZK under different real-world workloads.

## Requirement

### For AVX-MSM
- Ubuntu 22.04.4
- gcc version 11.4.0 
- cmake 3.22.1
- support AVX512-IFMA instruction sets

### For AVX-ZK
```
$ sudo apt-get install build-essential cmake git libgmp3-dev libprocps3-dev python-markdown libboost-all-dev libssl-dev
  ```
### For jsnark

- JDK 8 (Higher versions are also expected to work)
- Junit 4
- BouncyCastle library

## Build instructions

## AVX-MSM

### Building

 Target the `SimdMSM` library.

```shell
$ cd  AVX-MSM
$ cd demo/381/
$ make lib
```
If you encounter the error `../../../preset/x64-pbc-bls12-381.sh: not found`, try the following two commands: 
```shell
$ chmod +x ../../preset/x64-pbc-bls12-381.sh
$ sed -i 's/\r$//' ../../preset/x64-pbc-bls12-381.sh
```
### Using
Run AVX-MSM. The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pip_ifma.c`.

```shell
$ mkdir build
$ make ifma
$ ./build/test_pip_ifma
```

Run AVX-pair-MSM. The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pair_ifma.c`.

```shell
$ make pair_ifma
$ ./build/test_pair_ifma
```

Run AVX-MSM(muti-threads). The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pip_threads.c`.

```shell
$ make ifma
$ ./build/test_pip_threads
```

Run AVX-pair-MSM(muti-threads). The benchmark's data size `num` and window size `WSIZE` can be modified in the file `/test/test_pair_threads.c`.

```shell
$ make pair_ifma
$ ./build/test_pair_threads
```

You can also use the Python script to perform batch benching.
```shell
$ mkdir build
$ python bench.py
```
### Output example
The output structure of AVX-MSM, AVX-pair-MSM, AVX-MSM (multi-threads), and AVX-pair-MSM (multi-threads) is generally similar. Here, I'll use AVX-MSM as an example to describe its output structure.

The three macros `WNUM`, `WMBITS`, and `NBENCHS` in the test file represent the multi-scalar multiplication scale, window size, and number of benchmark iterations, respectively.

``` c
Pippenger_old=0.790256  // the execution time of the original Pippenger 
Pippenger_ifma=0.325606 // the execution time of our AVX-MSM (in seconds)
YES  // the computation result is correct
```

The output of `bench.py` is as follows: the first column represents the multi-scalar multiplication scale, followed by the window size, the execution time of the original Pippenger, the execution time of our AVX-MSM, and the speedup between the two.
```
[15, 6, 0.057, 0.019, 3.0]
[15, 7, 0.059, 0.019, 3.1052631578947367]
[15, 8, 0.058, 0.018, 3.2222222222222228]
[15, 9, 0.034, 0.014, 2.428571428571429]
[15, 10, 0.037, 0.017, 2.176470588235294]
[15, 11, 0.038, 0.02, 1.9]
[15, 12, 0.056, 0.027, 2.074074074074074]
[15, 13, 0.05, 0.034, 1.4705882352941175]
Best: [15, 9, 0.034, 0.014, 2.428571428571429] // this is the best window size
```

## AVX-ZK
### Building

Generate static link library `libmsm.a`.

```shell
$ cd AVX-MSM/demo/381
$ make msm
```

Cmake and create the Makefile:

```shell
$ cd AVX-ZK
$ mkdir build && cd build && cmake ..
```

Copy the `libmsm.a` and `librelic_s.a`.to libsnark/build/depends/libff/libff.
```shell
$ cp ../../AVX-MSM/demo/381/build/libmsm.a ../../AVX-MSM/demo/381/target/lib/librelic_s.a ./libsnark/build/depends/libff/libff
```

Then, to compile the library, run this within the `build` directory:

```shell
$ make
```

### Using

Run the profiling of AVX-ZK.

```shell
$ make profile_r1cs_gg_ppzksnark
$ ./libsnark/profile_r1cs_gg_ppzksnark 1000  100 bytes
```

### Running and Testing AVX-ZK by JsnarkCircuitBuilder

### Building
This part is similar to AVX-ZK.
```shell
$ cd jsnark/libsnark
$ mkdir build && cd build && cmake ..
```
Copy the `libmsm.a` and `librelic_s.a`.to libsnark/build/depends/libff/libff. Then build.
```shell
$ cp ../../../AVX-MSM/demo/381/build/libmsm.a ../../../AVX-MSM/demo/381/target/lib/librelic_s.a ./libsnark/build/depends/libff/libff
$ make
```

To compile the JsnarkCircuitBuilder project via command line, from the jsnark directory:

```shell
$ cd jsnark
$ cd JsnarkCircuitBuilder
$ mkdir -p bin
$ javac -d bin -cp /usr/share/java/junit4.jar:bcprov-jdk15on-159.jar  $(find ./src/* | grep ".java$")
```

### Using

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

