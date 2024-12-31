# AVX-ZK

This source code is an efficient implemetation of MSM and zkSNARK using AVX-512IFMA.

## Requirement

- Ubuntu 22.04.4
- gcc version 11.4.0 
- cmake 3.22.1
- support AVX512-IFMA instruction sets

## Build

### AVX-MSM

 Target the `relic `library.

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

