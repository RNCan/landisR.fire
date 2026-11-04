# Base Fire parameter compatibility layer

This development layer mirrors the parameter structures used by the current
LANDIS-II Base Fire source code while keeping them internal until strict
equivalence has been demonstrated.

## Reference structures

`FireRegion.cs` defines, by region:

- Name
- MapCode
- MeanSize
- MinSize
- MaxSize
- IgnitionProbability
- FireSpreadAge
- FuelCurve
- WindCurve
- Index

`InputParameterParser.cs` additionally reads:

- Timestep
- Species_CSV_File
- InitialFireRegionsMap
- optional DynamicFireRegionTable
- FuelCurveTable
- WindCurveTable
- FireDamageTable
- MapNames
- LogFile
- SummaryLogFile

## Compatibility choices

The R layer intentionally does not add a MinSize <= MeanSize <= MaxSize check,
because that relationship is not validated by the reference FireRegion/parser.

Region names and map codes are unique. Map codes are stored as unsigned-short
compatible values from 0 through 65535.

Fuel and wind curves are kept separate. Missing curve rows retain the C# object
default of zeros in the prepared cache.

FuelCurve uses strict increasing checks for non-negative thresholds after
severity 1. WindCurve uses the reference implementation's weaker rule that
rejects decreases but permits equality.

FireDamageTable requires increasing maximum-age proportions in (0, 1] and
strictly increasing severity-minus-tolerance values above -4 and no greater
than 4.

## Performance representation

`.basefire_prepare_parameter_cache()` converts validated data into:

- contiguous region vectors;
- 1-based internal region indices;
- retained 0-based Base Fire indices;
- a direct integer `map_code_to_region` lookup of length 65536;
- 5-column integer fuel and wind matrices;
- compact damage vectors;
- a named species fire-tolerance vector.

No data-frame join is required in later cell-to-cell fire kernels.

## Parser scope

The internal parser targets the current section names in the Base Fire source:
`DynamicFireRegionTable`, `FuelCurveTable`, `WindCurveTable` and
`FireDamageTable`. Legacy section aliases can be added later if required by
real historical parameter files.

The parser remains internal until strict Base Fire equivalence is demonstrated.
