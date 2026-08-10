# Dye issue investigation

This folder collects copies of the notebooks relevant to the passive-dye / GPU
kernel investigation. The original files remain in their existing locations.

## Contents

- `notebooks/current/amazon_river_debugging.ipynb`: current full debugging run.
- `notebooks/model/amazon_river_sim.ipynb`: original Amazon model notebook.
- `notebooks/dye_domainerror_minimal_reproducer.ipynb`: small, data-free
  CPU/GPU reduction ladder for the reported dye kernel exception.
- `notebooks/trial4-archive/`: archived Trial 4 variants 4.1 through 4.8.
- `tests/dye_domainerror_reproducer.jl`: standalone version suitable for
  attaching to an Oceananigans issue.
- `environment/`: the repository's Julia `Project.toml` and `Manifest.toml`.

Large JLD2 outputs are intentionally excluded. A minimal reproducer should not
require downloaded forcing data or existing simulation output.

## Next step

Open `notebooks/dye_domainerror_minimal_reproducer.ipynb` in a fresh Julia
kernel and run the four cases in order. The current host stacktrace reports the
exception during the `NaNChecker` synchronization, which may be later than the
kernel that actually failed. The reproducer enables launch blocking and
synchronizes after each timestep to localize the first failing case.

The standalone script can run one case at a time by setting `DYE_REPRO_CASE`
to `cpu_redi_immersed`, `gpu_no_redi_immersed`, `gpu_redi_regular`, or
`gpu_redi_immersed`. With no setting, it runs the full ladder.
