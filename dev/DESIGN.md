# landisR.fire design

## Dependency rule

`landisR.fire` may use `landisR`; `landisR` must never import or depend on `landisR.fire`.
The fire package is an external process model connected through a generic fire-event contract.

## Core decomposition

A fire model is decomposed into five replaceable modules:

1. ignition
2. target size / stopping rule
3. spatial spread
4. cell-level severity
5. mortality / forest impact

Post-fire ecological response such as serotiny remains the responsibility of landisR through `landis_fire()`.

## Base Fire strategy

The first scientific implementation will be a strict compatibility profile that reproduces LANDIS-II Base Fire.
It must use the same generic component interfaces as every later model. No Base Fire assumption should be hard-coded in the orchestration layer.

## Extension requirements

The architecture must permit, without modifying the core engine:

- empirical/custom fire-size frequency distributions;
- alternative continuous fire-size distributions;
- ignition and burn probability depending on forest composition, biomass, age, fuel, climate, FWI or user-defined predictors;
- cell-level severity depending on composition or other covariates;
- alternative spread algorithms;
- custom mortality relationships;
- completely user-defined component functions.

## Standard fire event

Every model should ultimately produce the same event representation: year, ignition cell, burned cells, spread order/source, cell-level severity, target area, realized area, fire region and metadata.
This standard object is the boundary between fire generation and application to the forest state.

## Development sequence

1. freeze the generic contracts;
2. implement strict Base Fire FireRegion + parameter parsing;
3. implement ignition;
4. implement target size;
5. implement spread and wind;
6. implement severity;
7. implement mortality;
8. connect fire_event to landisR::landis_fire();
9. validate event-by-event against LANDIS-II Base Fire;
10. add extended/custom modules only after strict compatibility is demonstrated.
