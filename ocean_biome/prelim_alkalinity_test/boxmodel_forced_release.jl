# second box test: add alkalinity as a timed forcing instead of an initial change
# both carbonate copies start the same and only Alk2 gets the release

using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField
using CSV

# short run so the release is resolved before, during, and after it happens
const TIMESTEP = 5seconds
const STOP_TIME = 2hours
const OUTPUT_INTERVAL = 5minutes
const RELEASE_START = 20minutes
const RELEASE_DURATION = 1hour

# background conditions, kept the same as the first box test for comparison
const TEMP = 25.0
const SALINITY = 35.0
const BG_DIC = 2100.0
const BG_ALK = 2300.0

# Alk2 should finish 100 mmol m^-3 above Alk1
const TARGET_ALK_ADDITION = 100.0
const ALK_RELEASE_RATE = TARGET_ALK_ADDITION / RELEASE_DURATION

const OUTPUT_DIRECTORY = joinpath(@__DIR__, "boxmodel_forced_release_outputs")
const JLD2_FILE = joinpath(OUTPUT_DIRECTORY, "box_forced_release.jld2")
const CSV_FILE = joinpath(OUTPUT_DIRECTORY, "box_forced_release.csv")
mkpath(OUTPUT_DIRECTORY)

# use the same biological starting point as the original box example
const INITIAL_BIOLOGY = (
    NO₃ = 10.0,
    NH₄ = 0.1,
    P = 0.1,
    Z = 0.01,
)

# keep the original seasonal light function even though this test is only two hours
const MODEL_YEAR = 365days
surface_PAR(t) = 60 * (1 - cos((t + 15days) * 2pi / MODEL_YEAR)) *
                 (1 / (1 + 0.2 * exp(-((mod(t, MODEL_YEAR) - 200days) / 50days)^2))) + 2

const BOX_DEPTH = -10.0
PAR_at_depth(t) = surface_PAR(t) * exp(0.2 * BOX_DEPTH)

# BoxModel forcing functions receive the model clock and its fields
# return a concentration tendency only during the release window
function alkalinity_release(clock, fields)
    t = clock.time
    during_release = RELEASE_START <= t < RELEASE_START + RELEASE_DURATION
    return during_release ? ALK_RELEASE_RATE : 0.0
end

# one box holds both carbonate copies so they experience identical biology
grid = BoxModelGrid()
clock = Clock(time = zero(grid))
PAR = FunctionField{Center, Center, Center}(PAR_at_depth, grid; clock)

biogeochemistry = LOBSTER(
    grid;
    inorganic_carbon = CarbonateSystem(2),
    open_bottom = false,
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR),
)

# only the second alkalinity tracer gets the timed source
model = BoxModel(;
    biogeochemistry,
    clock,
    forcing = (; Alk2 = alkalinity_release),
)

# control and treatment start from exactly the same carbonate state
set!(
    model;
    NO₃ = INITIAL_BIOLOGY.NO₃,
    NH₄ = INITIAL_BIOLOGY.NH₄,
    P = INITIAL_BIOLOGY.P,
    Z = INITIAL_BIOLOGY.Z,
    DIC1 = BG_DIC,
    Alk1 = BG_ALK,
    DIC2 = BG_DIC,
    Alk2 = BG_ALK,
)

simulation = Simulation(model; Δt = TIMESTEP, stop_time = STOP_TIME)
simulation.output_writers[:fields] = JLD2Writer(
    model,
    model.fields;
    filename = JLD2_FILE,
    schedule = TimeInterval(OUTPUT_INTERVAL),
    overwrite_existing = true,
)

progress(sim) = @info(
    "Forced-release progress",
    model_time = prettytime(time(sim)),
    alk2_minus_alk1 = model.fields.Alk2[1, 1, 1] - model.fields.Alk1[1, 1, 1],
)
simulation.callbacks[:progress] = Callback(progress, TimeInterval(20minutes))

@info(
    "Starting forced alkalinity box test",
    release_start = prettytime(RELEASE_START),
    release_duration = prettytime(RELEASE_DURATION),
    target_addition = TARGET_ALK_ADDITION,
)
run!(simulation)

# load both carbonate histories from the saved model output
DIC1_series = FieldTimeSeries(JLD2_FILE, "DIC1")
Alk1_series = FieldTimeSeries(JLD2_FILE, "Alk1")
DIC2_series = FieldTimeSeries(JLD2_FILE, "DIC2")
Alk2_series = FieldTimeSeries(JLD2_FILE, "Alk2")
times = DIC1_series.times

DIC1 = Float64.(DIC1_series[1, 1, 1, :])
Alk1 = Float64.(Alk1_series[1, 1, 1, :])
DIC2 = Float64.(DIC2_series[1, 1, 1, :])
Alk2 = Float64.(Alk2_series[1, 1, 1, :])

# calculate pH and pCO2 for the control and treatment at every saved time
carbon_chemistry = CarbonChemistry()
pH1 = similar(DIC1)
pH2 = similar(DIC2)
pCO2_1 = similar(DIC1)
pCO2_2 = similar(DIC2)

for n in eachindex(times)
    pH1[n] = carbon_chemistry(;
        DIC = DIC1[n], Alk = Alk1[n], T = TEMP, S = SALINITY,
        output = Val(:pHᶠ),
    )
    pH2[n] = carbon_chemistry(;
        DIC = DIC2[n], Alk = Alk2[n], T = TEMP, S = SALINITY,
        output = Val(:pHᶠ),
    )
    pCO2_1[n] = carbon_chemistry(;
        DIC = DIC1[n], Alk = Alk1[n], T = TEMP, S = SALINITY,
        output = Val(:pCO₂),
    )
    pCO2_2[n] = carbon_chemistry(;
        DIC = DIC2[n], Alk = Alk2[n], T = TEMP, S = SALINITY,
        output = Val(:pCO₂),
    )
end

# save a small comparison table that is easy to inspect and plot
CSV.write(
    CSV_FILE,
    (;
        time_minutes = times ./ minute,
        control_DIC = DIC1,
        control_Alk = Alk1,
        control_pH = pH1,
        control_pCO2 = pCO2_1,
        treatment_DIC = DIC2,
        treatment_Alk = Alk2,
        treatment_pH = pH2,
        treatment_pCO2 = pCO2_2,
        Alk_difference = Alk2 .- Alk1,
    ),
)

# basic bookkeeping checks for the release itself
alk_difference = Alk2 .- Alk1
before_release = times .< RELEASE_START
maximum_before_release = maximum(abs.(alk_difference[before_release]); init = 0.0)
final_addition = alk_difference[end]

isapprox(maximum_before_release, 0.0; atol = 1e-6) ||
    error("Control and treatment differ before the release")
isapprox(final_addition, TARGET_ALK_ADDITION; atol = 0.1) ||
    error("Final alkalinity addition was $final_addition instead of $TARGET_ALK_ADDITION")
all(isfinite, vcat(DIC1, Alk1, DIC2, Alk2, pH1, pH2, pCO2_1, pCO2_2)) ||
    error("A carbonate output contains NaN or Inf")

@info(
    "Forced alkalinity box test complete",
    maximum_difference_before_release = maximum_before_release,
    final_alkalinity_addition = final_addition,
    final_pH_change = pH2[end] - pH1[end],
    final_pCO2_change = pCO2_2[end] - pCO2_1[end],
    output_directory = OUTPUT_DIRECTORY,
)
