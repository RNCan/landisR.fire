# Base Fire target fire size compatibility

This layer reproduces `FireEvent.ComputeSize(IFireRegion ecoregion)`.

## Reference algorithm

For the current fire region:

1. if `MeanSize <= 0`, return 0 immediately;
2. set the exponential distribution rate to `1 / MeanSize`;
3. draw one continuous fire size;
4. if the draw is strictly below `MinSize`, use `MinSize`;
5. else if the draw is strictly above `MaxSize`, use `MaxSize`;
6. otherwise keep the draw;
7. divide by LANDIS-II `CellArea`;
8. cast to integer.

Because all valid fire sizes are non-negative, the C# integer cast is
equivalent to `floor()` after division by cell area.

## Important compatibility details

No minimum of one cell is imposed. A positive fire-size value smaller than one
cell area can therefore produce a target size of zero cells.

The Base Fire parser does not enforce `MinSize <= MeanSize <= MaxSize`.
`ComputeSize()` is therefore reproduced with its literal branch order rather
than silently sorting or correcting those parameters.

The lower and upper comparisons are strict (`<` and `>`). A draw exactly equal
to either boundary passes through the unbounded branch.

## RNG abstraction

The strict Base Fire size component accepts an `exponential_draw(lambda)`
function. The default R implementation uses `stats::rexp(1, rate=lambda)`.

This reproduces the statistical distribution and Base Fire control flow, but it
does not claim bit-for-bit equivalence with the LANDIS-II core random-number
generator. Deterministic tests inject known exponential draws so branching,
bounding, cell conversion and draw count can be verified independently.

A ModelCore-compatible RNG adapter can be added later if exact cross-engine RNG
stream equivalence is required.

## Prepared representation

The size component prepares:

- `mean_size`
- `min_size`
- `max_size`
- region indices
- cell area

into compact vectors/scalars. No regional join or parameter lookup is required
when a fire requests a target size.

## Future modular alternatives

This module is deliberately behind the generic `size` component contract.
Base Fire exponential size generation is therefore only one possible model.
The same orchestration can later use empirical frequencies, lognormal, gamma,
Weibull or fully custom target-size generators.
