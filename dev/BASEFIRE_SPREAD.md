# Base Fire spread, neighbors and wind compatibility

This layer reproduces the spatial control flow in `FireEvent.Spread()` and
`FireEvent.GetNeighbors()` while keeping Base Fire behind the generic
`spread` and `burn_probability` component contracts.

## Wind draws

At the start of each event Base Fire draws:

1. `windDirection = (int)(GenerateUniform() * 8)`
2. `windSpeed = GenerateUniform()`

`GetNeighbors()` clamps directions above 7 to 7.

For the eight first-ring directions N, NE, E, SE, S, SW, W and NW, the
reference wind probabilities are:

- primary and adjacent: `((4-speed)/8) * (1+speed)`
- perpendicular: `((4-speed)/8)`
- opposite and adjacent: `((4-speed)/8) * (1-speed)`

The probability vector is rotated using the literal C# index arithmetic.

Each geometrically existing first-ring site consumes one `NextDouble()` draw,
even if that site is inactive. A null/off-landscape neighbor consumes no draw.

A ninth candidate is two cells away in the wind direction. It consumes one
`GenerateUniform()` only when that geometric site exists, and is selected when
the draw is strictly below wind speed.

## Why Base Fire has its own site topology cache

The generic `fire_grid_topology()` intentionally contains active cells only.
That is ideal for most fast fire algorithms, but strict Base Fire RNG behavior
also depends on the existence of inactive sites.

The internal Base Fire topology therefore distinguishes:

- 0: no geometric site / null neighbor;
- positive full-site index with inactive flag;
- positive full-site index mapped to an active fire index.

Inactive candidates remain in the neighbor list through the shuffle and are
discarded only afterward, matching the reference code.

## Spread queue

The event runner uses an integer FIFO queue and logical membership vectors
rather than repeated list searches. This preserves `Queue<Site>` and
`Contains()` behavior while keeping lookup cost constant.

For each dequeued site the strict order is:

1. update current fire region and lazily compute its target size when needed;
2. assign the site to the event;
3. calculate severity;
4. only when severity > 0:
   - increment region and event fire-size counts;
   - update last-fire year;
   - apply the damage callback;
   - get wind-selected neighbors;
   - shuffle candidates;
   - skip inactive, already-event and already-queued sites;
   - apply the TSF / FireSpreadAge burn-probability gate;
   - enqueue successful candidates.

The event loop retains the original stopping conditions based on the current
region target and the maximum target encountered among fire regions.

## Separation of spread and burn probability

The Base Fire wind filter is the `spread` component.

The post-wind probability:

`1 - exp(-time_since_last_fire / FireSpreadAge)`

is the `burn_probability` component.

This separation is behaviorally equivalent to Base Fire but allows future
models to combine alternative spatial propagation algorithms with
composition-, fuel- or climate-sensitive conditional burning.

## RNG adapters

`GenerateUniform`, `NextDouble`, exponential size draws and shuffle are kept
behind separate injectable functions. The default R implementation does not
claim bit-for-bit equality with LANDIS-II ModelCore RNG or shuffle. Deterministic
tests verify control flow and draw consumption independently. A ModelCore RNG
adapter can be added later if exact random-stream equivalence is required.
