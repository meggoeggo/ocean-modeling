# reusable OAE pieces for the existing NumericalEarth Amazon model
# keep this file small: paired carbonate tracers plus one local Alk2 forcing

using OceanBioME
using Oceananigans
using Oceananigans.Grids: node
using Oceananigans.Units

struct AmazonOAEParameters{FT}
    longitude :: FT
    latitude  :: FT
    radius    :: FT
    depth     :: FT
    start_time :: FT
    duration   :: FT
    target_addition :: FT
end

function AmazonOAEParameters(;
    longitude = -49.5,
    latitude = 0.16,
    radius = 20kilometers,
    depth = 20meters,
    start_time = 1day,
    duration = 1day,
    target_addition = 100.0,
)
    values = promote(longitude, latitude, radius, depth,
                     start_time, duration, target_addition)
    return AmazonOAEParameters(values...)
end

const EARTH_RADIUS = 6.371e6
const DEG_TO_RAD = pi / 180

# approximate distance around the release center on the latitude-longitude grid
@inline function release_distance_squared(longitude, latitude, parameters)
    dy = (latitude - parameters.latitude) * EARTH_RADIUS * DEG_TO_RAD
    dx = (longitude - parameters.longitude) * EARTH_RADIUS * DEG_TO_RAD *
         cos(parameters.latitude * DEG_TO_RAD)
    return dx^2 + dy^2
end

# source tendency for treatment alkalinity in mmol m^-3 s^-1
@inline function amazon_alkalinity_release(i, j, k, grid, clock, model_fields, parameters)
    longitude, latitude, z = node(i, j, k, grid, Center(), Center(), Center())
    t = clock.time

    active_time = (parameters.start_time <= t) &
                  (t < parameters.start_time + parameters.duration)
    active_area = release_distance_squared(longitude, latitude, parameters) <= parameters.radius^2
    active_depth = (-parameters.depth <= z) & (z <= 0)

    rate = parameters.target_addition / parameters.duration
    return ifelse(active_time & active_area & active_depth, rate, zero(rate))
end

# create the LOBSTER model, paired carbonate tracers, and treatment forcing
function build_amazon_oae(grid; parameters = AmazonOAEParameters())
    biogeochemistry = LOBSTER(
        grid;
        inorganic_carbon = CarbonateSystem(2),
        open_bottom = false,
    )

    alkalinity_forcing = Forcing(
        amazon_alkalinity_release;
        discrete_form = true,
        parameters,
    )

    return (
        biogeochemistry = biogeochemistry,
        forcing = (; Alk2 = alkalinity_forcing),
        parameters = parameters,
    )
end

# the extra river closure uses a NamedTuple, so Oceananigans requires one
# diffusivity entry for every tracer after LOBSTER is added
function oae_river_diffusivities(river_kz)
    return (
        T = river_kz,
        S = river_kz,
        dye = river_kz,
        NO₃ = river_kz,
        NH₄ = river_kz,
        P = river_kz,
        Z = river_kz,
        DOM = river_kz,
        sPOM = river_kz,
        bPOM = river_kz,
        DIC1 = river_kz,
        DIC2 = river_kz,
        Alk1 = river_kz,
        Alk2 = river_kz,
        e = 0.0,
    )
end

# call this after ECCO has initialized T and S so those fields are not overwritten
function initialize_amazon_oae!(model;
    DIC = 2100.0,
    alkalinity = 2300.0,
    nitrate = 10.0,
    ammonium = 0.1,
    phytoplankton = 0.1,
    zooplankton = 0.01,
)
    set!(
        model;
        NO₃ = nitrate,
        NH₄ = ammonium,
        P = phytoplankton,
        Z = zooplankton,
        DOM = 0.0,
        sPOM = 0.0,
        bPOM = 0.0,
        DIC1 = DIC,
        Alk1 = alkalinity,
        DIC2 = DIC,
        Alk2 = alkalinity,
    )
    return nothing
end

function check_amazon_oae(model)
    tracer_names = keys(model.tracers)
    required = (:DIC1, :DIC2, :Alk1, :Alk2)
    all(name -> name in tracer_names, required) ||
        error("Amazon OAE model is missing one or more paired carbonate tracers")
    return true
end

function amazon_oae_outputs(model)
    tracers = model.tracers
    return (
        DIC1 = tracers.DIC1,
        DIC2 = tracers.DIC2,
        Alk1 = tracers.Alk1,
        Alk2 = tracers.Alk2,
        NO₃ = tracers.NO₃,
        NH₄ = tracers.NH₄,
        P = tracers.P,
        Z = tracers.Z,
        DOM = tracers.DOM,
        sPOM = tracers.sPOM,
        bPOM = tracers.bPOM,
    )
end
