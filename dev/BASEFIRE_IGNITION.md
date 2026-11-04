# Base Fire ignition compatibility

The first dynamic Base Fire kernel reproduces the initiation logic in
`FireEvent.Initiate()` and the time-since-fire probability in
`FireEvent.ComputeFireInitSpreadProb()`.

## Reference logic

For each active site, in landscape iteration order:

1. annual region ignition probability is multiplied directly by the fire
   extension timestep;
2. one uniform draw must be <= that value;
3. a second uniform draw must be <=
   `1 - exp(-time_since_last_fire / FireSpreadAge)`;
4. site severity is evaluated only after both probability gates pass;
5. severity must be > 0;
6. only then is a fire event initiated and spread immediately.

The extension does not clamp `IgnitionProbability * timestep` to one.

## RNG order

The R compatibility kernel deliberately evaluates one cell at a time.
The second uniform draw is not generated when the first gate fails. Severity is
also not evaluated until both probability gates pass.

This preserves the short-circuit structure of the C# code. A fully vectorized
landscape ignition draw would consume a different RNG stream and is therefore
not used for strict Base Fire compatibility.

## Initial TimeOfLastFire

`SiteVars.Initialize()` sets each active site's initial time of last fire to:

`StartTime - maximum cohort age at the site`

The ignition preparation layer reproduces this behavior when `last_fire_year`
is missing and `max_cohort_age` plus `context$start_time` are available.

## FireSpreadAge = 0

No additional guard is introduced. The literal reference formula is retained.
Therefore:

- positive time since fire with age 0 gives probability 1;
- zero time since fire with age 0 gives NaN.

A `uniform <= NaN` comparison fails, which is the effective C# behavior.

## Performance

Region indices, annual ignition probabilities, spread ages and last-fire years
are prepared once into compact vectors in active-cell order.

The future Base Fire run loop will remain sequential at the event-initiation
level because a successful event immediately spreads and updates
`TimeOfLastFire`, affecting later cells in the same landscape traversal.
Expensive joins or state extraction are still avoided inside this loop.
