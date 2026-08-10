# Adding OAE to the Amazon model

This version adds only the pieces already verified in the box experiment:

- LOBSTER
- two carbonate copies (`DIC1/Alk1` and `DIC2/Alk2`)
- identical control and treatment initialization
- one timed, localized forcing on `Alk2`

It intentionally does not add gas exchange yet. First verify that the regional
model transports the alkalinity perturbation without negative or non-finite
tracers.

## Changes to the Amazon run script

Load OceanBioME and include the component near the other imports:

```julia
using OceanBioME
include(joinpath(@__DIR__, "..", "oae", "amazon_oae.jl"))
```

After `grid` is constructed:

```julia
oae_parameters = AmazonOAEParameters(
    longitude = long_river,
    latitude = lat_river,
    radius = 20kilometers,
    depth = 20meters,
    start_time = 1day,
    duration = 1day,
    target_addition = 100.0,
)
oae = build_amazon_oae(grid; parameters = oae_parameters)
```

The custom river closure must now list every tracer. Replace its current `κ`
NamedTuple with:

```julia
river_mixing = VerticalScalarDiffusivity(
    VerticallyImplicitTimeDiscretization();
    κ = oae_river_diffusivities(river_mouth_kz),
)
```

Keep the existing dye sponge and add the alkalinity forcing:

```julia
model_forcing = merge(
    (; dye = dye_sponge),
    oae.forcing,
)
```

Add the biogeochemistry to the existing constructor:

```julia
ocean = ocean_simulation(
    grid;
    # existing arguments stay here
    tracers = (:T, :S, :dye),
    biogeochemistry = oae.biogeochemistry,
    forcing = model_forcing,
)
```

After this existing line initializes ECCO temperature and salinity:

```julia
set!(ocean.model, ecco_set)
```

initialize only the new biological and carbonate fields:

```julia
initialize_amazon_oae!(ocean.model)
check_amazon_oae(ocean.model)
```

Save the new fields with a 3D writer:

```julia
ocean.output_writers[:oae] = JLD2Writer(
    ocean.model,
    amazon_oae_outputs(ocean.model);
    filename = "amazon_oae_3d.jld2",
    schedule = TimeInterval(6hours),
    overwrite_existing = true,
)
```

## Before a scientific run

Confirm the release location, radius, depth, timing, and amount with your
mentor. `target_addition = 100` means that an unmixed selected cell receives a
100 mmol m^-3 concentration increase. It is not yet a total-moles or tonnes
release. A mass-based release must divide the total amount by the wet volume of
the selected cells on the immersed Amazon grid.
