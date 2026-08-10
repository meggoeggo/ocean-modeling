# check the alkalinity test results and make one summary figure
# this only reads saved files, so it does not rerun either model

using Oceananigans, Oceananigans.Units
using CairoMakie
using CSV

const HERE = @__DIR__
const ORIGINAL_BOX_FILE = normpath(joinpath(HERE, "..", "old_examples", "box.jld2"))

# look in both places since some of the older outputs were moved
const RESULT_DIRECTORIES = (
    joinpath(HERE, "boxmodel_edited_outputs"),
    HERE,
)

function find_result(filename)
    for directory in RESULT_DIRECTORIES
        path = joinpath(directory, filename)
        isfile(path) && return path
    end
    error("Could not find $filename in: $(join(RESULT_DIRECTORIES, ", "))")
end

# names for the files this script reads and creates
const CARBONATE_CSV = find_result("alkalinity_control_vs_treatment.csv")
const ANALYSIS_DIRECTORY = dirname(CARBONATE_CSV)
const REPORT_FILE = joinpath(ANALYSIS_DIRECTORY, "validation_report.txt")
const FIGURE_FILE = joinpath(ANALYSIS_DIRECTORY, "box_tracers_and_pco2.png")

# (1) broad ranges for now -- can modify later
const CONTEXT_BOUNDS = Dict(
    "DIC" => (1500.0, 2600.0),       # mmol m^-3
    "Alk" => (1800.0, 2800.0),       # mmol m^-3
    "pH" => (7.0, 8.5),              # free pH scale in this csv
    "pCO2" => (100.0, 2000.0),       # uatm
)

# (2) flag a value if it changes too much between saved points (also prelim)
const MAX_STEP_CHANGE = Dict(
    "DIC" => 100.0,                  # mmol m^-3 per saved interval
    "Alk" => 100.0,
    "pH" => 0.2,
    "pCO2" => 500.0,                 # uatm per saved interval
)

const EXPECTED_ALKALINITY_ADDITION = 100.0 # mmol m^-3
const ILLUSTRATIVE_ATMOSPHERIC_PCO2 = 420.0 # uatm, reference only

# store each check so it can be printed in one report later
function record!(results, level, check, detail)
    push!(results, (level = level, check = check, detail = detail))
end

function validate_carbonate_csv(path)
    # load each column using its csv header name
    data = CSV.File(path; types = Float64)
    isempty(data) && error("CSV is empty: $path")
    results = NamedTuple[]

    time_days = data.time_days
    control_DIC = data.control_DIC
    control_Alk = data.control_Alk
    control_pH = data.control_pH
    control_pCO2 = data.control_pCO2
    treatment_DIC = data.treatment_DIC
    treatment_Alk = data.treatment_Alk
    treatment_pH = data.treatment_pH
    treatment_pCO2 = data.treatment_pCO2

    # first check for broken numbers and incorrectly ordered output times
    all(name -> all(isfinite, getproperty(data, name)), propertynames(data)) ?
        record!(results, "PASS", "Finite values", "No NaN or Inf values") :
        record!(results, "FAIL", "Finite values", "At least one NaN or Inf value")

    all(diff(time_days) .> 0) ?
        record!(results, "PASS", "Time ordering", "Saved times increase strictly") :
        record!(results, "FAIL", "Time ordering", "Times repeat or decrease")

    # basic physical checks, not regional validation yet
    # for dic, alk, ph, pco2 -- can add the lobster values as well
    hard_checks = (
        ("control DIC", control_DIC, x -> x > 0, "must be positive"),
        ("treatment DIC", treatment_DIC, x -> x > 0, "must be positive"),
        ("control alkalinity", control_Alk, x -> x > 0, "must be positive in this seawater test"),
        ("treatment alkalinity", treatment_Alk, x -> x > 0, "must be positive in this seawater test"),
        ("control pH", control_pH, x -> 0 <= x <= 14, "must lie between 0 and 14"),
        ("treatment pH", treatment_pH, x -> 0 <= x <= 14, "must lie between 0 and 14"),
        ("control pCO2", control_pCO2, x -> x >= 0, "must be nonnegative"),
        ("treatment pCO2", treatment_pCO2, x -> x >= 0, "must be nonnegative"),
    )

    for (name, values, predicate, requirement) in hard_checks
        all(predicate, values) ?
            record!(results, "PASS", "Physical check: $name", requirement) :
            record!(results, "FAIL", "Physical check: $name", requirement)
    end

    # compare every variable against the broad ranges above (again, DIC/alk/ph/pco2)
    contextual = (
        ("control DIC", "DIC", control_DIC),
        ("treatment DIC", "DIC", treatment_DIC),
        ("control alkalinity", "Alk", control_Alk),
        ("treatment alkalinity", "Alk", treatment_Alk),
        ("control pH", "pH", control_pH),
        ("treatment pH", "pH", treatment_pH),
        ("control pCO2", "pCO2", control_pCO2),
        ("treatment pCO2", "pCO2", treatment_pCO2),
    )

    for (name, variable, values) in contextual
        lower, upper = CONTEXT_BOUNDS[variable]
        observed_min, observed_max = extrema(values)
        if observed_min >= lower && observed_max <= upper
            record!(results, "PASS", "Provisional range: $name",
                    "observed $observed_min to $observed_max; screen $lower to $upper")
        else
            record!(results, "WARN", "Provisional range: $name",
                    "observed $observed_min to $observed_max; screen $lower to $upper")
        end

        # look for sudden jumps between two saved ten-day points
        largest_step = maximum(abs.(diff(values)); init = 0.0)
        threshold = MAX_STEP_CHANGE[variable]
        level = largest_step <= threshold ? "PASS" : "WARN"
        record!(results, level, "Abrupt change: $name",
                "largest saved-step change $largest_step; provisional limit $threshold")
    end

    # DIC should match, while treatment Alk should stay 100 units higher
    dic_difference = treatment_DIC .- control_DIC
    alk_difference = treatment_Alk .- control_Alk
    maximum_DIC_mismatch = maximum(abs.(dic_difference))
    maximum_alk_offset_error = maximum(abs.(alk_difference .- EXPECTED_ALKALINITY_ADDITION))

    record!(results,
            maximum_DIC_mismatch <= 1e-6 ? "PASS" : "WARN",
            "Control/treatment DIC consistency",
            "maximum absolute DIC difference = $maximum_DIC_mismatch mmol m^-3")
    record!(results,
            maximum_alk_offset_error <= 1e-3 ? "PASS" : "WARN",
            "Alkalinity-offset consistency",
            "maximum deviation from +$(EXPECTED_ALKALINITY_ADDITION) = $maximum_alk_offset_error mmol m^-3")

    # save all pass/warn/fail results to a text file
    open(REPORT_FILE, "w") do io
        println(io, "Preliminary alkalinity experiment validation")
        println(io, "Source: $path")
        println(io, "Rows: $(length(time_days))")
        println(io, "Important: contextual bounds are provisional screening values.")
        println(io)
        for result in results
            println(io, "[$(result.level)] $(result.check): $(result.detail)")
        end
    end

    return (;
        time_days,
        control_pH,
        treatment_pH,
        control_pCO2,
        treatment_pCO2,
        results,
    )
end

validation = validate_carbonate_csv(CARBONATE_CSV)

# print a short summary in the terminal too
failures = count(result -> result.level == "FAIL", validation.results)
warnings = count(result -> result.level == "WARN", validation.results)
@info "Validation complete" failures warnings REPORT_FILE

# tracer names have to match the names stored in the original jld2 file
const BOX_TRACERS = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM)
const BOX_TITLES = ("NO3", "NH4", "P", "Z", "DOM", "sPOM")
const BOX_LABELS = (
    "Nitrate (mmol N m^-3)",
    "Ammonium (mmol N m^-3)",
    "Phytoplankton (mmol N m^-3)",
    "Zooplankton (mmol N m^-3)",
    "DOM (mmol N m^-3)",
    "Small POM (mmol N m^-3)",
)

# pull the original box tracer histories and convert time to years
box_series = NamedTuple{BOX_TRACERS}(
    FieldTimeSeries(ORIGINAL_BOX_FILE, string(tracer)) for tracer in BOX_TRACERS
)
box_times = first(values(box_series)).times ./ (365days)

# use the same two-column layout for all eight panels
fig = Figure(size = (1250, 1500), fontsize = 18)

# make six of the original LOBSTER tracer plots
for (index, tracer) in enumerate(BOX_TRACERS)
    row = div(index - 1, 2) + 1
    col = mod(index - 1, 2) + 1
    ax = Axis(
        fig[row, col];
        title = BOX_TITLES[index],
        xlabel = "Time (years)",
        ylabel = BOX_LABELS[index],
    )
    values_at_box = Float64.(box_series[tracer][1, 1, 1, :])
    lines!(ax, box_times, values_at_box; linewidth = 2.3)
end

# both carbonate plots use the saved csv time converted from days to years
carbonate_years = validation.time_days ./ 365

# replace the bPOM panel with the control and treatment pH
ax_ph = Axis(
    fig[4, 1];
    title = "Alkalinity treatment: pH response",
    xlabel = "Time (years)",
    ylabel = "pH (n/a)",
)
lines!(ax_ph, carbonate_years, validation.control_pH;
       linewidth = 2.5, label = "Control")
lines!(ax_ph, carbonate_years, validation.treatment_pH;
       linewidth = 2.5, label = "+100 mmol m^-3 Alk")
axislegend(ax_ph; position = :rt)

# show how much the treatment changes pH at the start and end
initial_delta_ph = validation.treatment_pH[1] - validation.control_pH[1]
final_delta_ph = validation.treatment_pH[end] - validation.control_pH[end]
text!(
    ax_ph,
    0.03,
    0.48;
    text = "Treatment - control: +$(round(initial_delta_ph; digits=3)) initially, " *
           "+$(round(final_delta_ph; digits=3)) finally",
    space = :relative,
    align = (:left, :center),
)


# pco2 comparison levels
ax_pco2 = Axis(
    fig[4, 2];
    title = "Alkalinity treatment: pCO2 response",
    xlabel = "Time (years)",
    ylabel = "pCO2 (uatm)",
)

# plot both cases against the same time axis
lines!(ax_pco2, carbonate_years, validation.control_pCO2;
       linewidth = 2.5, label = "Control")
lines!(ax_pco2, carbonate_years, validation.treatment_pCO2;
       linewidth = 2.5, label = "+100 mmol m^-3 Alk")
hlines!(ax_pco2, [ILLUSTRATIVE_ATMOSPHERIC_PCO2];
        linewidth = 1.8, linestyle = :dash, label = "Illustrative atmosphere")
axislegend(ax_pco2; position = :rt)

# calculate the treatment effect at the start and end and add it to the panel
initial_delta_pco2 = validation.treatment_pCO2[1] - validation.control_pCO2[1]
final_delta_pco2 = validation.treatment_pCO2[end] - validation.control_pCO2[end]
text!(
    ax_pco2,
    0.03,
    0.48;
    text = "Treatment - control: $(round(initial_delta_pco2; digits=1)) uatm initially, " *
           "$(round(final_delta_pco2; digits=1)) uatm finally",
    space = :relative,
    align = (:left, :center),
)

# save the finished figure beside the csv and report
save(FIGURE_FILE, fig)
@info "Saved validation figure" FIGURE_FILE
