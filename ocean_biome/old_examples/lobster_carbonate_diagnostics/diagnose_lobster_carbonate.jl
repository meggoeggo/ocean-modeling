# Diagnose biological divergence in LOBSTER + CarbonateSystem box experiments.
#
# Four cases isolate reproducibility, alkalinity sensitivity, and timestep
# sensitivity without changing box_model.jl or boxmodel_edited.jl.

using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField

const MODEL_YEAR = 365days
const STOP_TIME = MODEL_YEAR
const OUTPUT_INTERVAL = 1day
const AMBIENT_DIC = 2100.0
const AMBIENT_ALKALINITY = 2300.0
const ALKALINITY_ADDITION = 100.0
const NOMINAL_DEPTH = -10.0

const OUTPUT_DIRECTORY = joinpath(@__DIR__, "lobster_carbonate_diagnostics")
mkpath(OUTPUT_DIRECTORY)

const INITIAL_BIOLOGY = (
    NO₃ = 2.0,
    NH₄ = 0.05,
    P = 0.10,
    Z = 0.01,
    sPOM = 0.0,
    bPOM = 0.0,
    DOM = 0.0,
)

const TRACERS = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :DIC, :Alk)
const BIOLOGICAL_TRACERS = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM)

surface_PAR(t) = 60 * (1 - cos((t + 15days) * 2π / MODEL_YEAR)) *
                 (1 / (1 + 0.2 * exp(-((mod(t, MODEL_YEAR) - 200days) / 50days)^2))) + 2
PAR_at_depth(t) = surface_PAR(t) * exp(0.2 * NOMINAL_DEPTH)

function run_case(name; alkalinity, timestep)
    output_file = joinpath(OUTPUT_DIRECTORY, name * ".jld2")

    if isfile(output_file)
        @info "Reusing existing diagnostic output" name output_file
        return output_file
    end

    grid = BoxModelGrid()
    clock = Clock(time = zero(grid))
    PAR = FunctionField{Center, Center, Center}(PAR_at_depth, grid; clock)

    biogeochemistry = LOBSTER(
        grid;
        inorganic_carbon = CarbonateSystem(),
        open_bottom = false,
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR),
        # Disabled because the installed modifier expects model.tracers, while
        # BoxModel exposes model.fields.
        scale_negatives = false,
    )

    model = BoxModel(; biogeochemistry, clock)
    set!(model; INITIAL_BIOLOGY..., DIC = AMBIENT_DIC, Alk = alkalinity)

    simulation = Simulation(model; Δt = timestep, stop_time = STOP_TIME)
    simulation.output_writers[:fields] = JLD2Writer(
        model,
        model.fields;
        filename = output_file,
        schedule = TimeInterval(OUTPUT_INTERVAL),
        overwrite_existing = true,
    )

    @info "Running diagnostic case" name alkalinity timestep
    run!(simulation)
    return output_file
end

function load_case(file)
    series = NamedTuple{TRACERS}(FieldTimeSeries(file, string(tracer)) for tracer in TRACERS)
    times = first(Base.values(series)).times
    tracer_values = NamedTuple{TRACERS}(
        Float64.(series[tracer][1, 1, 1, :]) for tracer in TRACERS
    )
    return (; times, tracer_values...)
end

function max_difference(case_a, case_b, tracers)
    return maximum(
        maximum(abs.(case_a[tracer] .- case_b[tracer])) for tracer in tracers
    )
end

function differences_by_tracer(case_a, case_b)
    return NamedTuple{TRACERS}(
        maximum(abs.(case_a[tracer] .- case_b[tracer])) for tracer in TRACERS
    )
end

function first_difference_day(case_a, case_b, tracers; tolerance = 1e-10)
    for n in eachindex(case_a.times)
        difference = maximum(
            abs(case_a[tracer][n] - case_b[tracer][n]) for tracer in tracers
        )
        difference > tolerance && return case_a.times[n] / day
    end
    return missing
end

function write_timeseries_comparison(control, treatment, treatment_small_dt)
    csv_file = joinpath(OUTPUT_DIRECTORY, "diagnostic_differences.csv")

    open(csv_file, "w") do io
        header = ["time_days"]
        for tracer in TRACERS
            append!(header, ["delta_$(tracer)_alk", "delta_$(tracer)_small_dt"])
        end
        println(io, join(header, ','))

        for n in eachindex(control.times)
            row = Any[control.times[n] / day]
            for tracer in TRACERS
                push!(row, treatment[tracer][n] - control[tracer][n])
                push!(row, treatment_small_dt[tracer][n] - treatment[tracer][n])
            end
            println(io, join(row, ','))
        end
    end

    return csv_file
end


control_a_file = run_case(
    "control_a_dt5min";
    alkalinity = AMBIENT_ALKALINITY,
    timestep = 5minutes,
)
control_b_file = run_case(
    "control_b_dt5min";
    alkalinity = AMBIENT_ALKALINITY,
    timestep = 5minutes,
)
treatment_file = run_case(
    "treatment_dt5min";
    alkalinity = AMBIENT_ALKALINITY + ALKALINITY_ADDITION,
    timestep = 5minutes,
)
treatment_small_dt_file = run_case(
    "treatment_dt1min";
    alkalinity = AMBIENT_ALKALINITY + ALKALINITY_ADDITION,
    timestep = 1minute,
)

control_a = load_case(control_a_file)
control_b = load_case(control_b_file)
treatment = load_case(treatment_file)
treatment_small_dt = load_case(treatment_small_dt_file)

duplicate_difference = max_difference(control_a, control_b, TRACERS)
duplicate_biology_difference = max_difference(control_a, control_b, BIOLOGICAL_TRACERS)
alkalinity_biology_difference = max_difference(control_a, treatment, BIOLOGICAL_TRACERS)
timestep_biology_difference = max_difference(treatment, treatment_small_dt, BIOLOGICAL_TRACERS)
duplicate_differences_by_tracer = differences_by_tracer(control_a, control_b)
alkalinity_differences_by_tracer = differences_by_tracer(control_a, treatment)
timestep_differences_by_tracer = differences_by_tracer(treatment, treatment_small_dt)

first_alkalinity_divergence = first_difference_day(
    control_a,
    treatment,
    BIOLOGICAL_TRACERS;
    tolerance = 1e-10,
)
first_timestep_divergence = first_difference_day(
    treatment,
    treatment_small_dt,
    BIOLOGICAL_TRACERS;
    tolerance = 1e-10,
)

csv_file = write_timeseries_comparison(control_a, treatment, treatment_small_dt)

@info(
    "LOBSTER carbonate diagnostics complete",
    duplicate_max_abs_difference = duplicate_difference,
    duplicate_max_abs_biology_difference = duplicate_biology_difference,
    alkalinity_max_abs_biology_difference = alkalinity_biology_difference,
    timestep_max_abs_biology_difference = timestep_biology_difference,
    first_alkalinity_biology_divergence_day = first_alkalinity_divergence,
    first_timestep_biology_divergence_day = first_timestep_divergence,
    duplicate_differences_by_tracer,
    alkalinity_differences_by_tracer,
    timestep_differences_by_tracer,
    csv_file,
)
