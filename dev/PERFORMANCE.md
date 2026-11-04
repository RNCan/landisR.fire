# Performance strategy

`landisR.fire` is intended to remain highly modular while supporting large landscapes and many fire events.

## Rules

1. Extract forest state in bulk, never one `landis_state()` call per candidate cell.
2. Precompute regular-grid neighbor topology.
3. Use integer internal cell indices inside spread kernels.
4. Precompute or bulk-update composition/fuel susceptibility vectors when possible.
5. Avoid data.frame joins inside event-growth loops.
6. Keep event logs and user-facing tables outside the hottest kernels.
7. Preserve RNG order when strict Base Fire equivalence requires it.
8. Use vectorized R first and Rcpp only after profiling.

## Benchmark dimensions

Future benchmark scripts should report elapsed time separately for:

- preparation;
- ignition;
- target-size sampling;
- spatial spread;
- conditional burn probability;
- severity;
- mortality;
- conversion/application through the landisR interface.

Benchmarks should include both many-small-fire and few-large-fire scenarios because their bottlenecks may differ.
