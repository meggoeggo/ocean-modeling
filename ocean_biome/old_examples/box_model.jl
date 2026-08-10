
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField

# why?
const year = years = 365day

# perfectly mixed h20 parcel
grid = BoxModelGrid()

#set time to 0
clock = Clock(time = zero(grid))

# Prescribed time-dependent photosynthetically available radiation (PAR).
# PAR = photosynthetically available radiation --> light phytoplankton can use ofr photosynthesis
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * # creates annual light cycles
          (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2 #+2 --> keep par from being exactly 0

const z = -10 # Nominal depth of the box in meters.
PAR_func(t) = PAR⁰(t) * exp(0.2z) # control for light attenuatoin (light intensity dec as --> h20)
PAR = FunctionField{Center, Center, Center}(PAR_func, grid; clock)

# full box model compilarion
model = BoxModel(;
    biogeochemistry = LOBSTER( # LOBSTER: Lodyc-DAMTP Ocean Biogeochemical Simulation Tools for Ecosystem and Resources
        grid;
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR),
    ),
    clock,
)

set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01) # set init []s for four lobster tracers -- nitrate, ammonium, phytoplankton, zooplanton

# construct output file
output_file = joinpath(@__DIR__, "box.jld2")
simulation = Simulation(model; Δt = 5minutes, stop_time = 5years)
simulation.output_writers[:fields] = JLD2Writer(
    model,
    model.fields;
    filename = output_file,
    schedule = TimeInterval(10days),
    overwrite_existing = true,
)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"
simulation.callbacks[:progress] = Callback(prog, IterationInterval(1_000_000))

@info "Running the model..."
run!(simulation)

# Load all saved tracer time series.
times = FieldTimeSeries(output_file, "P").times
timeseries = NamedTuple{keys(model.fields)}(
    FieldTimeSeries(output_file, string(field))[1, 1, 1, :]
    for field in keys(model.fields)
)

# Plot the tracer time series. Makie layouts are one-based, so add 1 to the
# row and column calculated from the zero-based loop position.
using CairoMakie

fig = Figure(size = (1200, 1200), fontsize = 24)
axs = Axis[]

for (name, tracer) in pairs(timeseries)
    idx = length(axs)
    row = div(idx, 2) + 1
    col = mod(idx, 2) + 1
    ax = Axis(fig[row, col]; ylabel = string(name), xlabel = "Year", xticks = 0:5)
    push!(axs, ax)
    lines!(ax, times / year, tracer; linewidth = 3)
end

figure_file = joinpath(@__DIR__, "box_results.png")
save(figure_file, fig)
@info "Saved box-model figure" figure_file
