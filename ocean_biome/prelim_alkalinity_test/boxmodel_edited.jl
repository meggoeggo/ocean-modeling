
# this code attempts to implement alkalinity addition in a biogeochemical model
# 1- control- ambient DIC and total alkalinity
# 2- treatment- ambient DIC plus a known alkalinity addition
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField

# add consts
const MODEL_YR = 365days
const STOP_TIME = 5MODEL_YR
const TIMESTEP = 5minutes
const OUTPUT_INTERVAL = 10days

# backgroun tracer constants
const TEMP = 25.0
const SALINITY = 35.0
const BG_DIC = 2100.0 # identical initial DIC in both cases; changes during sim
const BG_ALK = 2300.0 # control alkalinity; test case adds 100 mmol m^-3 for 2400 alkalinty

# the actual added material
const ALK_ADDITION = 100.0 # mmol m⁻³

# outputs
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "boxmodel_edited_outputs")
mkpath(OUTPUT_DIRECTORY)

const INITIAL_BIOLOGY = ( # same as box model
    NO₃ = 10.0,
    NH₄ = 0.1,
    P = 0.1,
    Z = 0.01,
)

# par - kept same
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / MODEL_YR)) * # same as old model -- surface par
                 (1 / (1 + 0.2 * exp(-((mod(t, MODEL_YR) - 200days) / 50days)^2))) + 2

const z = -10.0 # m depth
PAR_func(t) = PAR⁰(t) * exp(0.2 * z) # par @ depth


# run case
function run_case(case_name, init_alk; stop_time = STOP_TIME, write_output = true)
    output_file = joinpath(OUTPUT_DIRECTORY, case_name * ".jld2")

    # define grid, time, PAR
    grid = BoxModelGrid()
    clock = Clock(time = zero(grid))
    PAR = FunctionField{Center, Center, Center}(PAR_func, grid; clock)

    # define LOBSTER setup
    biogeochemistry = LOBSTER(
        grid;
        inorganic_carbon = CarbonateSystem(), # explicitly integrate a carbonate system
        open_bottom = false, # make sure that the matter doesnt "fall" out of the box
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR),
    )

    # define boxmodel
    model = BoxModel(; biogeochemistry, clock)

    #set init params
    set!(
        model;
        INITIAL_BIOLOGY...,
        DIC = BG_DIC, # make actual model, integrate pre-set biogeochemistry with DIC and ALK
        Alk = init_alk,
    )

    simulation = Simulation(model; Δt = TIMESTEP, stop_time)

    if write_output
        simulation.output_writers[:fields] = JLD2Writer(
            model,
            model.fields;
            filename = output_file,
            schedule = TimeInterval(OUTPUT_INTERVAL),
            overwrite_existing = true,
        )
    end

    progress(sim) = @info(
        "Case progress",
        case = case_name,
        model_time = prettytime(time(sim)),
        wall_time = prettytime(simulation.run_wall_time),
    )
    simulation.callbacks[:progress] = Callback(progress, TimeInterval(MODEL_YR))

    @info "Starting carbonate case" case_name init_alk write_output
    run!(simulation)
    @info "Completed carbonate case" case_name

    return output_file
end

# calculate prelim carbonate diagnostic values
function carbonate_diagnostics(output_file)
    DIC = FieldTimeSeries(output_file, "DIC")
    Alk = FieldTimeSeries(output_file, "Alk")
    times = DIC.times
    # pull the values

    # convert to proper units
    carbon_chemistry = CarbonChemistry()
    pH = zeros(length(times))
    pCO₂ = zeros(length(times))

    for n in eachindex(times)
        dic = Float64(DIC[1, 1, 1, n])
        alk = Float64(Alk[1, 1, 1, n])

        # calculate ph at given timestep
        pH[n] = carbon_chemistry(;
            DIC = dic,
            Alk = alk,
            T = TEMP,
            S = SALINITY,
            output = Val(:pHᵗ),
        )
        # calc pC02 at given timestep
        pCO₂[n] = carbon_chemistry(;
            DIC = dic,
            Alk = alk,
            T = TEMP,
            S = SALINITY,
            output = Val(:pCO₂),
        )
    end

    return (; times, DIC, Alk, pH, pCO₂)
end

# save outputs as csv files
function write_comparison_csv(control, treatment)
    csv_file = joinpath(OUTPUT_DIRECTORY, "alkalinity_control_vs_treatment.csv")

    open(csv_file, "w") do io
        println(
            io,
            "time_days,control_DIC,control_Alk,control_pH,control_pCO2," *
            "treatment_DIC,treatment_Alk,treatment_pH,treatment_pCO2",
        )

        for n in eachindex(control.times)
            control_DIC = control.DIC[1, 1, 1, n]
            control_Alk = control.Alk[1, 1, 1, n]
            treatment_DIC = treatment.DIC[1, 1, 1, n]
            treatment_Alk = treatment.Alk[1, 1, 1, n]

            println(
                io,
                join(
                    (
                        control.times[n] / day,
                        control_DIC,
                        control_Alk,
                        control.pH[n],
                        control.pCO₂[n],
                        treatment_DIC,
                        treatment_Alk,
                        treatment.pH[n],
                        treatment.pCO₂[n],
                    ),',',
                ),
            )
        end
    end

    @info "Saved carbonate comparison" csv_file
    return csv_file
end

# create the final sims/files
control_file = run_case("carbonate_control", BG_ALK)
treatment_file = run_case(
    "alkalinity_treatment",
    BG_ALK + ALK_ADDITION,
)

control = carbonate_diagnostics(control_file)
treatment = carbonate_diagnostics(treatment_file)
write_comparison_csv(control, treatment)

@info(
    "Carbonate experiment complete",
    ALK_ADDITION = ALK_ADDITION,
    output_directory = OUTPUT_DIRECTORY,
)
