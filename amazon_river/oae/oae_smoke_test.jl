# tiny test of the exact NumericalEarth + river mixing + OAE interface

using NumericalEarth
using Oceananigans
using Oceananigans.TurbulenceClosures: VerticalScalarDiffusivity,
                                          VerticallyImplicitTimeDiscretization

include(joinpath(@__DIR__, "amazon_oae.jl"))

grid = RectilinearGrid(
    CPU();
    size = (8, 8, 4),
    halo = (5, 5, 4),
    x = (0, 1),
    y = (0, 1),
    z = (-40, 0),
)

# x and y are treated like degrees here only to test the forcing interface
parameters = AmazonOAEParameters(
    longitude = 0.5,
    latitude = 0.5,
    radius = 100kilometers,
    depth = 40meters,
    start_time = 0.0,
    duration = 100.0,
    target_addition = 1.0,
)
oae = build_amazon_oae(grid; parameters)

test_kz(x, y, z, t) = 1e-3
vertical_mixing = NumericalEarth.Oceans.default_ocean_closure()
river_mixing = VerticalScalarDiffusivity(
    VerticallyImplicitTimeDiscretization();
    κ = oae_river_diffusivities(test_kz),
)

ocean = ocean_simulation(
    grid;
    Δt = 1.0,
    closure = (vertical_mixing, river_mixing),
    tracers = (:T, :S, :dye),
    biogeochemistry = oae.biogeochemistry,
    forcing = oae.forcing,
    radiative_forcing = nothing,
    momentum_advection = WENOVectorInvariant(order = 5),
    tracer_advection = WENO(order = 5),
    verbose = false,
)

set!(ocean.model, T = 25.0, S = 35.0, dye = 0.0)
initialize_amazon_oae!(ocean.model)
check_amazon_oae(ocean.model)

ocean.stop_time = 10.0
run!(ocean)

alk_difference = maximum(ocean.model.tracers.Alk2) - maximum(ocean.model.tracers.Alk1)
isapprox(alk_difference, 0.1; atol = 1e-6) ||
    error("Expected an Alk difference of 0.1, got $alk_difference")

@info "Amazon OAE smoke test passed" alk_difference keys(ocean.model.tracers)
