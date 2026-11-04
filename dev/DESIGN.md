# landisR.fire design

## Dependency rule

`landisR.fire` may use `landisR`; `landisR` must never import or depend on `landisR.fire`.
The fire package is an external process model connected through a generic fire-event contract.

## Core decomposition

A fire model is decomposed into six replaceable modules:

1. ignition
2. target size / stopping rule
3. spatial spread / exposure
4. conditional burn probability
5. cell-level severity
6. mortality / forest impact

The explicit separation between spread and burn probability is intentional. Spread determines which candidate cells are reachable/exposed. Burn probability determines whether an exposed candidate actually burns. This allows, for example, Base Fire spatial propagation to be combined later with composition-sensitive flammability without rewriting the spread algorithm.

Post-fire ecological response such as serotiny remains the responsibility of landisR through `landis_fire()`.

## Prepare once / execute fast

Scientific modularity must not imply repeated expensive work inside cell-level spread loops.

The runtime architecture therefore has two phases:

1. configuration and preparation;
2. event execution using compact prepared structures.

`fire_prepare()` calls optional component-specific preparation functions once. These may create vectors, matrices, lookup tables or other reusable caches. Runtime components should consume those caches directly.

For regular raster landscapes, `fire_grid_topology()` precomputes:

- active-cell internal integer indices;
- mapping from internal index to external `cell_id`;
- 8 first-ring neighbors;
- 8 second-ring directional neighbors;
- compact row/column lookup geometry.

Missing or inactive neighbors are encoded as integer 0. Internal indices are used during spread because they are faster and smaller than repeated joins on arbitrary cell identifiers.

## Landscape-state access

The fire engine must never call `landis_state()` once per candidate cell. Required forest information should be extracted in bulk at the beginning of a fire timestep or other appropriate update interval and converted into a compact fire landscape view.

Potential cached predictors include:

- total biomass;
- species biomass/composition;
- cohort-derived age metrics;
- fire region index;
- time since fire;
- fuel variables;
- FWI/climate variables;
- topographic covariates.

A model that uses composition-dependent flammability should normally calculate or update its cell-level susceptibility vector in bulk, then index that vector during spread.

## Base Fire strategy

The first scientific implementation will be a strict compatibility profile that reproduces LANDIS-II Base Fire.
It must use the same generic component interfaces as every later model. No Base Fire assumption should be hard-coded in the orchestration layer.

Conceptually:

```r
fire_model_basefire()
```

will correspond to the strict combination of Base Fire ignition, size, spread, burn probability, severity and mortality components.

## Extension requirements

The architecture must permit, without modifying the core engine:

- empirical/custom fire-size frequency distributions;
- alternative continuous fire-size distributions;
- ignition depending on forest composition, biomass, age, fuel, climate, FWI or user-defined predictors;
- conditional burn probability depending on those same predictors;
- cell-level severity depending on composition or other covariates;
- alternative spread algorithms;
- custom mortality relationships;
- completely user-defined component functions.

## Standard fire event

Every model should ultimately produce the same event representation: year, ignition cell, burned cells, spread order/source, cell-level severity, target area, realized area, fire region and metadata.
This standard object is the boundary between fire generation and application to the forest state.

## Performance policy

Optimization follows measurement:

1. reproduce the intended biology;
2. validate against the reference implementation or controlled tests;
3. profile ignition, size, spread, burn probability, severity, mortality and landisR interface time separately;
4. optimize only demonstrated hotspots;
5. rerun equivalence tests after every optimization.

Vectorized R and precomputation are preferred first. Rcpp is reserved for kernels where profiling shows a meaningful benefit.

## Development sequence

1. freeze generic contracts including burn probability and preparation;
2. implement strict Base Fire FireRegion + parameter parsing;
3. implement ignition;
4. implement target size;
5. implement spread and wind;
6. implement Base Fire conditional burning rules;
7. implement severity;
8. implement mortality;
9. connect fire_event to landisR::landis_fire();
10. validate event-by-event against LANDIS-II Base Fire;
11. profile and optimize hotspots;
12. add extended/custom modules only after strict compatibility is demonstrated.


## Base Fire parameter layer

The Base Fire compatibility profile keeps its parameter representation separate
from the generic fire-model contracts.

The internal parameter layer stores:

- timestep;
- species fire tolerance;
- fire-region name and map code;
- mean, minimum and maximum event size;
- ignition probability;
- fire spread age;
- initial and optional dynamic fire-region maps;
- fuel severity curves;
- wind severity curves;
- fire damage classes;
- map/log output names.

Preparation converts these tables into compact numeric vectors, integer lookup
arrays and matrices. In particular, map codes use a direct 65,536-entry integer
lookup so event kernels can convert a raster map code to a 1-based internal
region index without joins or searches.

No additional ordering constraint is imposed on Base Fire MinSize, MeanSize and
MaxSize because the reference FireRegion/InputParameterParser does not impose
one. FuelCurve and WindCurve validation also intentionally remain distinct
because the reference C# setters use slightly different ordering rules.

These are compatibility details of the Base Fire profile, not generic
requirements of landisR.fire.

## Base Fire ignition kernel

Strict Base Fire compatibility uses sequential event initiation in active-cell
landscape order.

The ignition component separates two probability concepts already present in
the reference C#:

1. annual `IgnitionProbability * timestep`;
2. `ComputeFireInitSpreadProb`, based on time since last fire and
   `FireSpreadAge`.

The second random draw occurs only if the first test succeeds. Severity is
evaluated only if both random-probability tests succeed, and the site must have
severity > 0 before an event is constructed.

This short-circuit order is preserved even though a vectorized implementation
would be shorter, because pre-generating all random values would change the RNG
stream. Event initiation also cannot be safely precomputed for a whole
timestep: a successful fire spreads immediately and updates last-fire state
before later active cells are visited.

Performance is obtained instead by preparing region indices, annual ignition
probabilities, fire-spread ages and last-fire years into contiguous vectors.

## Base Fire target-size component

`FireEvent.ComputeSize()` is implemented as the first concrete target-size
component.

The compatibility kernel preserves:

- immediate zero when `MeanSize <= 0`;
- exponential rate `1 / MeanSize`;
- one distribution draw otherwise;
- strict `< MinSize` and `> MaxSize` bounding;
- literal branch behavior even for unordered Min/Mean/Max inputs;
- conversion from area to number of cells by division by `CellArea` followed
  by integer truncation;
- no forced minimum of one cell.

The exponential generator is injected behind a small function boundary. This
keeps the Base Fire statistical model independent of the package RNG
implementation and will let a ModelCore-compatible generator be substituted if
cross-engine random-stream equivalence is required.

This module is also the first demonstration of the intended extension design:
the generic fire engine depends on a `size` component, not on Base Fire's
exponential distribution.

## Base Fire spread and burn-probability components

Base Fire spatial behavior is decomposed without changing its control flow.

`spread` reproduces `GetNeighbors()`:

- 8-direction first-ring geometry;
- wind-direction rotation;
- wind-speed probability vector;
- separate `NextDouble()` draws for existing first-ring sites;
- optional ninth site two cells downwind.

`burn_probability` reproduces the later
`GenerateUniform() <= ComputeFireInitSpreadProb()` gate.

A Base Fire-specific full-site topology is prepared in addition to the generic
active-cell topology because inactive sites still consume `NextDouble()` and
participate in the shuffled candidate list in the reference C#.

The event runner uses an integer FIFO queue plus logical `in_queue` and
`event_occupied` vectors. This is semantically equivalent to the C#
`Queue<Site>` / `Contains()` checks but avoids repeated linear searches.

The runner currently receives severity and damage as callbacks. This keeps the
spatial layer testable before severity and mortality are implemented while also
preserving the intended modular architecture.
