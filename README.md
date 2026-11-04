# landisR.fire

Development package for modular fire models that interact with landisR.

The first scientific target is strict reproduction of LANDIS-II Base Fire. Base Fire will be implemented as one configuration of replaceable modules for ignition, target size, spread, conditional burn probability, severity and mortality.

The architecture is intentionally open so users can later substitute empirical fire-size distributions, composition-sensitive burning/severity, alternative spread algorithms or completely custom functions without modifying landisR.

## Performance design

Scientific flexibility is separated from runtime execution. Components may define a one-time preparation function, and `fire_prepare()` stores reusable numeric caches before events are simulated. Regular-grid neighbor topology is precomputed with internal integer indices through `fire_grid_topology()` so cell-to-cell spread does not need repeated joins or spatial searches.

The intended pattern is:

```r
model <- fire_model(
  ignition = ...,
  size = ...,
  spread = ...,
  burn_probability = ...,
  severity = ...,
  mortality = ...
)

prepared <- fire_prepare(model, landscape)
```

Custom components can remain pure R functions, while built-in components may later use vectorized R or Rcpp kernels where profiling demonstrates a real bottleneck.

## Base Fire compatibility work

The first Base Fire-specific layer is now kept internal while equivalence is being established. It represents and validates the native Base Fire region parameters, dynamic region schedule, fuel curves, wind curves, damage classes and species fire tolerance. It also prepares compact numeric lookup structures for later ignition, spread, severity and mortality kernels.

Base Fire-specific functions are deliberately not exported yet. They will become public only after event-level behavior has been validated against LANDIS-II Base Fire.
