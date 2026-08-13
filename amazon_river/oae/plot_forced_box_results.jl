# plot every evolving tracer and carbonate diagnostic from the forced box experiment

using Oceananigans
using Oceananigans.Units
using CSV
using CairoMakie

# find the output files made by boxmodel_forced_release.jl
const FORCED_OUTPUT_DIRECTORY = normpath(joinpath(
    @__DIR__,
    "..",
    "..",
    "ocean_biome",
    "prelim_alkalinity_test",
    "boxmodel_forced_release_outputs",
))

const JLD2_FILE = joinpath(FORCED_OUTPUT_DIRECTORY, "box_forced_release.jld2")
const CSV_FILE = joinpath(FORCED_OUTPUT_DIRECTORY, "box_forced_release.csv")
const FIGURE_FILE = joinpath(FORCED_OUTPUT_DIRECTORY, "box_forced_release_all_results.png")

# these match the release timing used in the model file
const RELEASE_START_MINUTES = 20.0
const RELEASE_END_MINUTES = 80.0

isfile(JLD2_FILE) || error("Model output was not found: $JLD2_FILE")
isfile(CSV_FILE) || error("Carbonate output was not found: $CSV_FILE")

# read a box-model field and turn its saved values into a simple vector
function read_box_field(field_name)
    field = FieldTimeSeries(JLD2_FILE, field_name)
    times_minutes = field.times ./ minute
    values = Float64.(field[1, 1, 1, :])
    return times_minutes, values
end

# show when the alkalinity forcing starts and stops on each plot
function mark_release!(axis)
    vlines!(
        axis,
        [RELEASE_START_MINUTES, RELEASE_END_MINUTES];
        color = :gray50,
        linestyle = :dash,
        linewidth = 1.2,
    )
end

# read the carbonate values calculated after the model run
carbonate = CSV.File(CSV_FILE; types = Float64)
carbonate_time = collect(carbonate.time_minutes)

# list the biological fields exactly as they are named in the JLD2 output
biology_fields = (
    ("NO₃", "NO3"),
    ("NH₄", "NH4"),
    ("P", "Phytoplankton"),
    ("Z", "Zooplankton"),
    ("DOM", "Dissolved organic matter"),
    ("sPOM", "Small particulate organic matter"),
    ("bPOM", "Large particulate organic matter"),
)

fig = Figure(size = (1500, 1600), fontsize = 17)
Label(
    fig[0, 1:3],
    "Timed alkalinity forcing: complete box-model output";
    fontsize = 25,
)

axes = Axis[]

# plot every biological tracer saved by LOBSTER
for (index, (field_name, title)) in enumerate(biology_fields)
    row = div(index - 1, 3) + 1
    column = mod(index - 1, 3) + 1
    axis = Axis(
        fig[row, column];
        title,
        xlabel = "Time (minutes)",
        ylabel = "Concentration (mmol m^-3)",
    )

    times, values = read_box_field(field_name)
    lines!(axis, times, values; linewidth = 2.5)
    mark_release!(axis)
    push!(axes, axis)
end

# compare control and treatment DIC
ax_dic = Axis(
    fig[3, 2];
    title = "DIC",
    xlabel = "Time (minutes)",
    ylabel = "DIC (mmol m^-3)",
)
lines!(ax_dic, carbonate_time, carbonate.control_DIC; linewidth = 2.5, label = "Control")
lines!(ax_dic, carbonate_time, carbonate.treatment_DIC; linewidth = 2.5, label = "Forced")
mark_release!(ax_dic)
axislegend(ax_dic; position = :rb)
push!(axes, ax_dic)

# compare control and treatment alkalinity
ax_alk = Axis(
    fig[3, 3];
    title = "Alkalinity",
    xlabel = "Time (minutes)",
    ylabel = "Alkalinity (mmol m^-3)",
)
lines!(ax_alk, carbonate_time, carbonate.control_Alk; linewidth = 2.5, label = "Control")
lines!(ax_alk, carbonate_time, carbonate.treatment_Alk; linewidth = 2.5, label = "Forced")
mark_release!(ax_alk)
axislegend(ax_alk; position = :rb)
push!(axes, ax_alk)

# isolate the alkalinity added by the forcing
ax_difference = Axis(
    fig[4, 1];
    title = "Added alkalinity",
    xlabel = "Time (minutes)",
    ylabel = "Forced - control (mmol m^-3)",
)
lines!(ax_difference, carbonate_time, carbonate.Alk_difference; linewidth = 2.5)
hlines!(ax_difference, [100.0]; color = :gray50, linestyle = :dot, linewidth = 1.2)
mark_release!(ax_difference)
push!(axes, ax_difference)

# plot the carbonate response to the alkalinity release
ax_ph = Axis(
    fig[4, 2];
    title = "pH response",
    xlabel = "Time (minutes)",
    ylabel = "pH",
)
lines!(ax_ph, carbonate_time, carbonate.control_pH; linewidth = 2.5, label = "Control")
lines!(ax_ph, carbonate_time, carbonate.treatment_pH; linewidth = 2.5, label = "Forced")
mark_release!(ax_ph)
axislegend(ax_ph; position = :rb)
push!(axes, ax_ph)

ax_pco2 = Axis(
    fig[4, 3];
    title = "pCO2 response",
    xlabel = "Time (minutes)",
    ylabel = "pCO2 (uatm)",
)
lines!(ax_pco2, carbonate_time, carbonate.control_pCO2; linewidth = 2.5, label = "Control")
lines!(ax_pco2, carbonate_time, carbonate.treatment_pCO2; linewidth = 2.5, label = "Forced")
mark_release!(ax_pco2)
axislegend(ax_pco2; position = :rt)
push!(axes, ax_pco2)

# use one shared time range so every panel lines up with the release period
linkxaxes!(axes...)
xlims!(first(axes), minimum(carbonate_time), maximum(carbonate_time))

save(FIGURE_FILE, fig)
@info "Saved complete forced-release figure" FIGURE_FILE

fig
