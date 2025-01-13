##  AVX-MSM

The main codes of `AVX-MSM` is in the `AVX-MSM/demo/381/` directory.

**No multi-threading:**
-  `pip_ifma` is MSM implementation using AVX512-IFMA in $\mathbb{G}_{1}$ group.
-  `pair_ifma` is MSM implementation using AVX512-IFMA in both ($\mathbb{G}_{1}$, $\mathbb{G}_{2}$) group.

**With multi-threading:**
-  `pip_threads` is MSM implementation using AVX512-IFMA in $\mathbb{G}_{1}$ group.
-  `pair_threads` is MSM implementation using AVX512-IFMA in both ($\mathbb{G}_{1}$, $\mathbb{G}_{2}$) group.

**Others:**
- `bench.py` is a python script for large-scale data testing.

