ENV["CUDA_LAUNCH_BLOCKING"] = "1"

using Pkg
using Oceananigans
using Oceananigans.Units
using CUDA
using Printf
using Oceananigans.TurbulenceClosures:
    IsopycnalSkewSymmetricDiffusivity,
    AdvectiveFormulation

CUDA.allowscalar(false)

const Nx = 8
const Ny = 8
const Nz = 8
const Lx = 100kilometers
const Ly = 100kilometers
const H  = 1000meters
const Δt = 30seconds
const number_of_steps = 10

"""Construct a tiny dye model with optional GM/Redi and immersed bathymetry."""
function build_dye_model(arch; use_redi, immersed)
    z = ExponentialDiscretization(Nz, -H, 0; scale=H/4)

    underlying_grid = RectilinearGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = (3, 3, 3),
        x = (0, Lx),
        y = (0, Ly),
        z,
        topology = (Bounded, Bounded, Bounded),
    )

    if immersed
        # A smooth seamount leaves shallow and deep wet columns without data files.
        bottom(x, y) = -900meters +
                       550meters * exp(-((x - Lx/2)^2 + (y - Ly/2)^2) / (20kilometers)^2)
        grid = ImmersedBoundaryGrid(
            underlying_grid,
            GridFittedBottom(bottom);
            active_cells_map = true,
        )
    else
        grid = underlying_grid
    end

    closure = if use_redi
        IsopycnalSkewSymmetricDiffusivity(
            κ_skew = 1e3,
            κ_symmetric = 1e3,
            skew_flux_formulation = AdvectiveFormulation(),
        )
    else
        nothing
    end

    model = HydrostaticFreeSurfaceModel(
        grid;
        buoyancy = BuoyancyTracer(),
        tracers = (:b, :dye),
        closure,
        momentum_advection = WENOVectorInvariant(order=5),
        tracer_advection = WENO(order=5),
        free_surface = SplitExplicitFreeSurface(grid; substeps=20),
    )

    # Stable stratification with a small horizontal slope for GM/Redi.
    bᵢ(x, y, z) = 1e-5 * z + 1e-7 * (x - Lx/2)

    # A strictly nonnegative near-surface passive dye patch.
    dyeᵢ(x, y, z) = exp(-((x - Lx/2)^2 + (y - Ly/2)^2) / (12kilometers)^2) *
                     exp(-((z + 40meters)^2) / (30meters)^2)

    set!(model, b=bᵢ, dye=dyeᵢ)
    return model
end

"""Run one case, synchronizing after every step to localize GPU failures."""
function run_case(name, arch; use_redi, immersed)
    @info "START" name arch use_redi immersed
    model = build_dye_model(arch; use_redi, immersed)

    initial_dye = Array(interior(model.tracers.dye))
    @assert minimum(initial_dye) >= 0
    @assert all(isfinite, initial_dye)

    minimum_dye = minimum(initial_dye)
    maximum_dye = maximum(initial_dye)

    for step in 1:number_of_steps
        time_step!(model, Δt)
        arch isa GPU && CUDA.synchronize()

        dye = Array(interior(model.tracers.dye))
        @assert all(isfinite, dye) "non-finite dye at step $step"
        minimum_dye = min(minimum_dye, minimum(dye))
        maximum_dye = max(maximum_dye, maximum(dye))
    end

    @info "PASS" name minimum_dye maximum_dye
    return (; name, minimum_dye, maximum_dye)
end

function report_environment()
    println("Julia version: ", VERSION)
    Pkg.status(["Oceananigans", "CUDA"])
    if CUDA.functional()
        println("CUDA.functional() = true")
    else
        println("CUDA.functional() = false")
    end
end

const available_cases = (
    cpu_redi_immersed = ("CPU + GM/Redi + immersed", CPU(), true, true),
    gpu_no_redi_immersed = ("GPU + no GM/Redi + immersed", GPU(), false, true),
    gpu_redi_regular = ("GPU + GM/Redi + regular", GPU(), true, false),
    gpu_redi_immersed = ("GPU + GM/Redi + immersed", GPU(), true, true),
)

function run_named_case(case_name::Symbol)
    haskey(available_cases, case_name) ||
        error("Unknown DYE_REPRO_CASE=$case_name. Choose one of $(keys(available_cases)).")

    name, arch, use_redi, immersed = available_cases[case_name]
    arch isa GPU && !CUDA.functional() && error("CUDA.functional() = false")
    return run_case(name, arch; use_redi, immersed)
end

function run_reproducer()
    report_environment()

    requested_case = Symbol(get(ENV, "DYE_REPRO_CASE", "all"))
    if requested_case != :all
        return [run_named_case(requested_case)]
    end

    results = []
    push!(results, run_case("CPU + GM/Redi + immersed", CPU(); use_redi=true, immersed=true))

    if !CUDA.functional()
        @warn "CUDA is unavailable; GPU cases skipped"
        return results
    end

    # Controls run first. The suspected case runs last because a device-side
    # exception may require restarting Julia before another GPU run.
    push!(results, run_case("GPU + no GM/Redi + immersed", GPU(); use_redi=false, immersed=true))
    push!(results, run_case("GPU + GM/Redi + regular", GPU(); use_redi=true, immersed=false))
    push!(results, run_case("GPU + GM/Redi + immersed", GPU(); use_redi=true, immersed=true))

    return results
end

try
    results = run_reproducer()
    println("\nAll requested cases completed.")
    foreach(println, results)
catch err
    println(stderr, "\nREPRODUCER FAILED")
    showerror(stderr, err, catch_backtrace())
    println(stderr)
    rethrow()
end
