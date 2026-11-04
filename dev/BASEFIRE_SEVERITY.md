# Base Fire severity compatibility

This layer reproduces `FireEvent.CalcSeverity()`.

## Fuel severity

For the site's current fire region, Base Fire calculates time since last fire:

`currentTime - TimeOfLastFire[site]`

It then evaluates severity thresholds 1 through 5 in order. A threshold applies
when it is not -1 and time since last fire is greater than or equal to it.

Later matching thresholds overwrite earlier values, so the final fuel severity
is the highest matching severity.

## Wind severity

Wind severity is evaluated only when the external LANDIS-II site variable
`Wind.TimeOfLastEvent` exists.

Time since last wind is:

`currentTime - TimeOfLastWind[site]`

Wind thresholds are also evaluated from severity 1 through 5, but the
comparison is reversed:

`timeSinceLastWind <= threshold`

Again, later matching severities overwrite earlier values.

## Final severity

The final site severity is:

`max(fuelSeverity, windSeverity)`

and therefore lies from 0 through 5 under valid Base Fire parameters.

## External wind state

In LANDIS-II, Base Fire does not create `TimeOfLastWind`; it requests the
registered site variable `Wind.TimeOfLastEvent`. The R compatibility layer
therefore treats wind history as optional external state.

The preparation layer accepts it from either:

- `fire_state$cells$last_wind_year`, or
- `context$time_of_last_wind`.

When neither exists, wind severity is skipped entirely even though wind curves
may be present.

## Performance

Region indices, fire history, optional wind history, fuel curves and wind
curves are prepared once. Per-cell severity then requires only direct vector and
matrix indexing.

The calculation has no random draws.
