# Base Fire cohort mortality compatibility

This layer reproduces the cohort filter used by
`FireEvent.IDisturbance.ReduceOrKillMarkedCohort()`.

## Severity 5

Fire severity 5 kills every cohort immediately.

The reference code takes this branch before species fire tolerance, longevity or
damage-table calculations. The R kernel preserves that short circuit.

## Severities 0 through 4

For lower severities Base Fire calculates:

`severity_tolerance_difference = siteSeverity - speciesFireTolerance`

and:

`age_as_percent = cohortAge / speciesLongevity`

The damage table is then scanned until a row satisfies both:

`severity_tolerance_difference >= SeverTolerDifference`

and:

`MaxAge >= age_as_percent`

The cohort is killed at the first matching row.

If no row matches, the cohort survives.

## Damage-table order

The current parser requires damage rows to have strictly increasing `MaxAge`
and strictly increasing `SeverTolerDifference`, and appends them to
`FireDamages` in that input order.

`PlugIn.Initialize()` passes `parameters.FireDamages` directly to
`FireEvent.Initialize()`.

The mortality callback then uses `foreach (IDamageTable damage in damages)`.
Therefore the executed behavior is a forward scan of the increasing table.

An old inline comment says "Read table backwards; most severe first", but the
current code does not reverse the list. The compatibility kernel follows the
executed code, not that stale comment.

## Float precision

The C# callback explicitly casts cohort age, species longevity and damage
`MaxAge` to `float` for the age comparison.

The R compatibility kernel converts those quantities to IEEE-754 single
precision before comparison so exact boundary behavior is preserved as closely
as possible.

## Removal fraction

A killed cohort returns `1.0`; a surviving cohort returns `0.0`.

There is no partial biomass reduction in this Base Fire callback even though
the Universal Cohort disturbance interface supports fractional reduction.

## Species longevity

Fire tolerance comes from the Base Fire species CSV/cache.

Longevity is not a Base Fire input parameter in LANDIS-II. It comes from the
core species dataset through `cohort.Species.Longevity`.

Accordingly, the R mortality preparation layer accepts longevity from external
model context through either:

- a named `context$species_longevity` vector, or
- `context$species` with `species_code` and `longevity`.

This keeps Base Fire dependent on landisR/core species state rather than
duplicating longevity inside fire parameters.
