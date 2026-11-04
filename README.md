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

## Base Fire ignition

The internal Base Fire compatibility profile now includes the event-initiation
kernel from `FireEvent.Initiate()`. It preserves the two-stage probability
gate, time-since-fire response, severity > 0 gate and short-circuit RNG order.

The kernel is intentionally evaluated one active cell at a time because a
successful Base Fire event spreads immediately and changes fire state before
later landscape cells are tested. Region-level values and last-fire years are
prepared into compact vectors first, so this sequential control flow does not
require repeated data extraction or joins.

Base Fire remains internal until full event-level equivalence is validated.

## Base Fire target size

The internal Base Fire profile now reproduces `FireEvent.ComputeSize()`.
`MeanSize` parameterizes an exponential distribution through
`lambda = 1 / MeanSize`; the generated area is bounded by `MinSize` and
`MaxSize`, divided by cell area, and truncated to an integer number of cells.

A non-positive `MeanSize` returns zero without consuming a size draw. Base Fire
also does not force positive generated areas to at least one raster cell.

The implementation remains behind the generic replaceable `size` component, so
future empirical or alternative fire-size distributions do not require changes
to ignition, spread, severity or mortality.

## Base Fire spatial spread

The internal Base Fire profile now includes the wind-driven `GetNeighbors()`
logic and the FIFO event-spread control flow.

The eight nearest geometric sites are filtered with the original directional
wind probabilities, and a ninth site two cells downwind may also be selected.
For strict compatibility, inactive geometric sites remain in the candidate list
through neighbor selection and shuffling, then are skipped before the
time-since-fire gate.

The original post-wind spread probability is implemented as the separate
`burn_probability` component. This means Base Fire now demonstrates four
concrete replaceable process modules: ignition, size, spread and conditional
burn probability.

Severity and mortality remain injected callbacks/placeholders until their own
compatibility kernels are implemented and validated.

## Base Fire severity

The internal Base Fire profile now reproduces `FireEvent.CalcSeverity()`.

Fuel severity increases with time since the last fire according to the region's
five `FuelCurve` thresholds. If an external wind-disturbance history is
available, wind severity is also calculated using the region's `WindCurve`
thresholds and time since the last wind event. Final severity is the maximum of
the fuel and wind severities.

Wind history remains optional because Base Fire itself does not create
`Wind.TimeOfLastEvent`; LANDIS-II obtains that site variable from another
extension when available.

Severity is still a replaceable process component, so future models can use
composition, biomass, fuel types, weather or other predictors without changing
Base Fire ignition, size or spatial spread.

## Base Fire cohort mortality

The internal Base Fire profile now reproduces the cohort-damage callback used
by `FireEvent`.

Severity 5 kills all cohorts. At lower severities, mortality depends on the
difference between site fire severity and species fire tolerance, together with
cohort age as a fraction of species longevity and the Base Fire damage table.

A matching cohort is removed completely (`1.0`); otherwise it survives
unchanged (`0.0`). Base Fire does not perform partial biomass reduction in this
callback.

Species fire tolerance remains part of the Base Fire parameter cache, while
species longevity is supplied by external landisR/core species context, matching
the division of responsibility in LANDIS-II.

All six fundamental Base Fire process components are now internally
implemented. Complete event-level assembly and cross-engine validation remain
to be done before Base Fire is exposed as a public profile.
