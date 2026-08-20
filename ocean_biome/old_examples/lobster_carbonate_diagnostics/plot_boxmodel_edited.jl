# Plot diagnostics from the existing carbonate control and alkalinity-treatment
# runs without rerunning either simulation.

using OceanBioME, Oceananigans, Oceananigans.Units
using CairoMakie

const TEMPERATURE = 25.0 # °C
const SALINITY = 35.0    # PSU

const OUTPUT_DIRECTORY = joinpath(@__DIR__, "boxmodel_edited_outputs")
const CONTROL_FILE = joinpath(OUTPUT_DIRECTORY, "carbonate_control.jld2")
const TREATMENT_FILE = joinpath(OUTPUT_DIRECTORY, "alkalinity_treatment.jld2")

for file in (CONTROL_FILE, TREATMENT_FILE)
    isfile(file) || error("Missing model output: $file")
end

"""Load one scalar tracer from a box-model JLD2 output."""
function load_tracer(file, tracer)
    series = FieldTimeSeries(file, tracer)
    return series.times, Float64.(series[1, 1, 1, :])
end

times, control_DIC = load_tracer(CONTROL_FILE, "DIC")
_, treatment_DIC = load_tracer(TREATMENT_FILE, "DIC")
_, control_Alk = load_tracer(CONTROL_FILE, "Alk")
_, treatment_Alk = load_tracer(TREATMENT_FILE, "Alk")
_, control_P = load_tracer(CONTROL_FILE, "P")
_, treatment_P = load_tracer(TREATMENT_FILE, "P")
_, control_NO₃ = load_tracer(CONTROL_FILE, "NO₃")
_, treatment_NO₃ = load_tracer(TREATMENT_FILE, "NO₃")

carbon_chemistry = CarbonChemistry()

function diagnose_carbonate(DIC, Alk)
    pH = similar(DIC)
    pCO₂ = similar(DIC)

    for n in eachindex(DIC)
        pH[n] = carbon_chemistry(;
            DIC = DIC[n], Alk = Alk[n], T = TEMPERATURE, S = SALINITY,
            output = Val(:pHᵗ),
        )
        pCO₂[n] = carbon_chemistry(;
            DIC = DIC[n], Alk = Alk[n], T = TEMPERATURE, S = SALINITY,
            output = Val(:pCO₂),
        )
    end

    return pH, pCO₂
end

control_pH, control_pCO₂ = diagnose_carbonate(control_DIC, control_Alk)
treatment_pH, treatment_pCO₂ = diagnose_carbonate(treatment_DIC, treatment_Alk)
time_years = times ./ (365days)

fig = Figure(size = (1200, 1250), fontsize = 19)

function comparison_axis(position; title, ylabel)
    return Axis(position; title, xlabel = "Time (years)", ylabel)
end

ax_alk = comparison_axis(fig[1, 1]; title = "Total alkalinity", ylabel = "Alk (mmol m⁻³)")
lines!(ax_alk, time_years, control_Alk; linewidth = 2.5, label = "Control")
lines!(ax_alk, time_years, treatment_Alk; linewidth = 2.5, label = "+100 mmol m⁻³ Alk")
axislegend(ax_alk; position = :rb)

ax_dic = comparison_axis(fig[1, 2]; title = "Dissolved inorganic carbon", ylabel = "DIC (mmol m⁻³)")
lines!(ax_dic, time_years, control_DIC; linewidth = 2.5, label = "Control")
lines!(ax_dic, time_years, treatment_DIC; linewidth = 2.5, label = "Treatment")
axislegend(ax_dic; position = :rb)

ax_ph = comparison_axis(fig[2, 1]; title = "Carbonate response: pH", ylabel = "pH (free scale)")
lines!(ax_ph, time_years, control_pH; linewidth = 2.5, label = "Control")
lines!(ax_ph, time_years, treatment_pH; linewidth = 2.5, label = "Treatment")
axislegend(ax_ph; position = :rb)

ax_pco2 = comparison_axis(fig[2, 2]; title = "Carbonate response: pCO₂", ylabel = "pCO₂ (µatm)")
lines!(ax_pco2, time_years, control_pCO₂; linewidth = 2.5, label = "Control")
lines!(ax_pco2, time_years, treatment_pCO₂; linewidth = 2.5, label = "Treatment")
axislegend(ax_pco2; position = :rt)

ax_carbon_delta = comparison_axis(
    fig[3, 1];
    title = "Treatment minus control",
    ylabel = "Concentration difference (mmol m⁻³)",
)
lines!(ax_carbon_delta, time_years, treatment_Alk .- control_Alk;
       linewidth = 2.5, label = "ΔAlk")
lines!(ax_carbon_delta, time_years, treatment_DIC .- control_DIC;
       linewidth = 2.5, label = "ΔDIC")
axislegend(ax_carbon_delta; position = :rc)

ax_biology_delta = comparison_axis(
    fig[3, 2];
    title = "Biological consistency check",
    ylabel = "Treatment − control (mmol N m⁻³)",
)
lines!(ax_biology_delta, time_years, treatment_P .- control_P;
       linewidth = 2.5, label = "ΔP")
lines!(ax_biology_delta, time_years, treatment_NO₃ .- control_NO₃;
       linewidth = 2.5, label = "ΔNO₃")
axislegend(ax_biology_delta; position = :rc)

figure_file = joinpath(OUTPUT_DIRECTORY, "alkalinity_control_vs_treatment.png")
save(figure_file, fig)

@info(
    "Saved carbonate diagnostic figure",
    figure_file,
    final_delta_alkalinity = treatment_Alk[end] - control_Alk[end],
    final_delta_DIC = treatment_DIC[end] - control_DIC[end],
    maximum_abs_delta_P = maximum(abs.(treatment_P .- control_P)),
    maximum_abs_delta_NO₃ = maximum(abs.(treatment_NO₃ .- control_NO₃)),
)
