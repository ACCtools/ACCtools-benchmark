# ACCtools benchmark

This repository was created to benchmark the SKYPE pipeline using the
supplementary data released with the Severus manuscript. The benchmark uses
the published cell-line data, structural-variant calls, and truth sets to
evaluate SKYPE alongside other structural-variant callers.

The supporting dataset is available on Zenodo:
[Supporting data for the Severus manuscript](https://zenodo.org/records/14541057).

`skype_bench.py` runs native SKYPE first for each selected cell line and saves
its successful graph-limit pair as `limit_combinations.json`. Every subsequent
VCF-mode case receives that file through stage 02's `--limit_combinations`
option. A VCF case that cannot build the graph with the native pair is marked
failed; stage 02 does not try a fallback pair.
