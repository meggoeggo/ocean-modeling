# Adapted from the official OceanBioME one-dimensional column example:
# https://oceanbiome.github.io/OceanBioME.jl/stable/generated/column/

using OceanBioME, Oceananigans, Printf
using CairoMakie
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Units

const year = years = 365days

# Idealized North Atlantic surface light and mixed-layer forcing.
@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) *
                            (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)
@inline fmld1(t) = H(t, 50days, year) *
                   (1 / (1 + exp(-(t - 100days) / 5days))) *
                   (1 / (1 + exp((t - 330days) / 25days)))
@inline MLD(t) = -(10 + 340 * (1 - fmld1(year - eps(year)) *
                               exp(-mod(t, year) / 25days) -
                               fmld1(mod(t, year))))
@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4
@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10

grid = RectilinearGrid(
    size = (1, 1, 50),
    extent = (20meters, 20meters, 200meters),
)

biogeochemistry = LOBSTER(
    grid;
    surface_PAR = PAR⁰,
    inorganic_carbon = CarbonateSystem(),
    scale_negatives = true,
)

CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()
clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = ConstantField(35.0)

model = NonhydrostaticModel(
    grid;
    clock,
    closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
    biogeochemistry,
    boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),),
    auxiliary_fields = (; T, S),
)

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05,
            DIC = 2239.8, Alk = 2409.0)

simulation = Simulation(model; Δt = 3minutes, stop_time = 100days)
progress(sim) = @printf(
    "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
    iteration(sim),
    prettytime(sim),
    prettytime(sim.Δt),
    prettytime(sim.run_wall_time),
)
simulation.callbacks[:progress] = Callback(progress, TimeInterval(10days))

output_prefix = joinpath(@__DIR__, "column")
simulation.output_writers[:profiles] = JLD2Writer(
    model,
    model.tracers;
    filename = output_prefix * ".jld2",
    schedule = TimeInterval(1day),
    overwrite_existing = true,
)

qCO₂ = BoundaryConditionOperation(model.tracers.DIC, :top, model)
simulation.output_writers[:carbon_flux] = JLD2Writer(
    model,
    (; qCO₂);
    indices = (:, :, grid.Nz),
    filename = output_prefix * "_carbon.jld2",
    schedule = TimeInterval(1day),
    overwrite_existing = true,
)

@info "Running one-dimensional LOBSTER column model"
run!(simulation)

# Load the saved output and diagnose sinking carbon export.
P = FieldTimeSeries(output_prefix * ".jld2", "P")
NO₃ = FieldTimeSeries(output_prefix * ".jld2", "NO₃")
Z = FieldTimeSeries(output_prefix * ".jld2", "Z")
sPOM = FieldTimeSeries(output_prefix * ".jld2", "sPOM")
bPOM = FieldTimeSeries(output_prefix * ".jld2", "bPOM")
DIC = FieldTimeSeries(output_prefix * ".jld2", "DIC")
Alk = FieldTimeSeries(output_prefix * ".jld2", "Alk")
air_sea_CO₂_flux = FieldTimeSeries(output_prefix * "_carbon.jld2", "qCO₂")

_, _, z = nodes(P)
times = P.times
carbon_export = zeros(length(times))
R_CN = model.biogeochemistry.underlying_biogeochemistry.plankton.carbon_ratio
export_k = grid.Nz - 20

for (n, t) in enumerate(times)
    clock.time = t
    sPOM_flux = sPOM[n][1, 1, export_k] *
                biogeochemical_drift_velocity(model.biogeochemistry, Val(:sPOM)).w[1, 1, export_k]
    bPOM_flux = bPOM[n][1, 1, export_k] *
                biogeochemical_drift_velocity(model.biogeochemistry, Val(:bPOM)).w[1, 1, export_k]
    carbon_export[n] = (sPOM_flux + bPOM_flux) * R_CN
end

# Plot the biological tracers and carbon fluxes.
fig = Figure(size = (1000, 1500), fontsize = 20)
axis_kwargs = (
    xlabel = "Time (days)",
    ylabel = "z (m)",
    limits = ((0, times[end] / days), (-150meters, 0)),
)

axP = Axis(fig[1, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(axP, times / days, z, interior(P, 1, 1, :, :)'; colormap = :batlow)
Colorbar(fig[1, 2], hmP)

axNO₃ = Axis(fig[2, 1]; title = "Nitrate concentration (mmol N / m³)", axis_kwargs...)
hmNO₃ = heatmap!(axNO₃, times / days, z, interior(NO₃, 1, 1, :, :)'; colormap = :batlow)
Colorbar(fig[2, 2], hmNO₃)

axZ = Axis(fig[3, 1]; title = "Zooplankton concentration (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(axZ, times / days, z, interior(Z, 1, 1, :, :)'; colormap = :batlow)
Colorbar(fig[3, 2], hmZ)

axD = Axis(fig[4, 1]; title = "Detritus concentration (mmol N / m³)", axis_kwargs...)
detritus = interior(sPOM, 1, 1, :, :)' .+ interior(bPOM, 1, 1, :, :)'
hmD = heatmap!(axD, times / days, z, detritus; colormap = :batlow)
Colorbar(fig[4, 2], hmD)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol
ax_flux = Axis(
    fig[5, 1];
    xlabel = "Time (days)",
    ylabel = "Flux (kgCO₂/m²/year)",
    title = "Air-sea CO₂ flux and sinking export",
    limits = ((0, times[end] / days), nothing),
)
lines!(ax_flux, times / days,
       interior(air_sea_CO₂_flux, 1, 1, 1, :) ./ 1e3 * CO₂_molar_mass * year * 10;
       linewidth = 3, label = "Air-sea flux ×10")
lines!(ax_flux, times / days, carbon_export / 1e3 * CO₂_molar_mass * year;
       linewidth = 3, label = "Sinking export")
Legend(fig[5, 2], ax_flux; framevisible = false)

figure_path = joinpath(@__DIR__, "column_results.png")
save(figure_path, fig)
@info "Saved column-model figure" figure_path
