# Complete internal Base Fire profile

This layer assembles the six independently validated Base Fire process kernels
into one strict internal profile.

## Six-component model

The profile is constructed with the generic `fire_model()` contract:

1. ignition
2. target size
3. spatial spread
4. conditional burn probability
5. severity
6. mortality

No Base Fire process is hard-coded into the generic `fire_model` constructor.

## Preparation

`.basefire_prepare_profile()` builds one shared Base Fire parameter cache and
passes it through the generic `fire_prepare()` mechanism.

All six process caches are prepared before timestep execution.

The complete geometric landscape, including inactive sites, is retained only in
the Base Fire spread cache because strict GetNeighbors() RNG behavior requires
it.

## Timestep orchestration

The internal `.basefire_run_timestep()` mirrors `PlugIn.Run()`:

1. reset current-timestep event ownership, severity and disturbed maps;
2. update the fire-region map if a dynamic region change is scheduled;
3. traverse active sites in topology/landscape order;
4. attempt Base Fire ignition at each site;
5. immediately spread every successful event;
6. update fire history and cohort state before evaluating later sites;
7. build event and summary diagnostics.

This sequencing is essential because earlier events can change
`TimeOfLastFire` before later active sites are evaluated.

## Cohort mortality

For the internal integration harness, cohorts are supplied as a data.frame with
at least:

- `cell_id`
- `species_code`
- `age`

The mortality component is called through the generic model contract for every
cohort at a disturbed site. Killed cohorts are removed immediately so a later
event in the same timestep sees the updated forest state.

This is an internal validation representation. The subsequent landisR handoff
will replace this table mutation with the package's generic fire application
interface.

## Three event concepts

LANDIS-II Base Fire contains several related but not identical notions of an
event.

The R result therefore keeps:

- `raw_events`: every successful FireEvent initiation;
- `events`: standard `landisR_fire_event` objects for events with at least one
  burned cell;
- `event_log`: rows that satisfy Base Fire's `FireEvent.Severity > 0` logging
  rule.

The summary `number_fires` follows the original plug-in and counts every
non-null FireEvent, while `total_sites_burned` and region totals follow the
event-log path.

This distinction preserves unusual legacy behavior without contaminating the
generic fire-event contract.

## Dynamic fire regions

Dynamic region updates are supported through an injected
`context$dynamic_region_map_reader(map_name, current_time, landscape)`.

The reader may return either:

- one map-code vector per landscape row, or
- a data.frame with `cell_id` and `fire_region_map_code`.

When a scheduled change occurs, all region-dependent prepared caches are
rebuilt once before the active-site loop.

## RNG control

The complete runner keeps separate injectable controls for:

- ignition uniforms;
- Base Fire `GenerateUniform`;
- Base Fire `NextDouble`;
- exponential fire-size draws;
- candidate shuffle;
- optional fixed wind direction and speed.

This makes deterministic event-level tests possible without claiming
bit-for-bit equivalence with LANDIS-II ModelCore RNG streams.

## Public API status

The complete profile remains internal. Public exposure should wait until:

1. full event-level behavior is validated;
2. the handoff to `landisR::landis_fire()` is implemented;
3. cross-engine validation against LANDIS-II Base Fire is satisfactory.
