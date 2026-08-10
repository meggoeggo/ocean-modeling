# Amazon River salinity diagnostics: Trial Set 4

These notebooks isolate the source of the negative-salinity behavior in the
Amazon River configuration. Each retained simulation notebook includes its
configuration and recorded output. Raw JLD2 fields and MP4 animations are not
versioned.

| Trial | Change tested | Recorded salinity result |
| --- | --- | --- |
| 4.1 | Spread river forcing over a 3 x 3 stencil | Day 12 minimum `-2.604316` |
| 4.2 | Add local upper-30-m vertical mixing | Day 12 minimum `-2.604316` |
| 4.3 | Locate and rank the actual river-forcing cells | Diagnostic only |
| 4.4 | Move the mixing region onto the measured forcing cells | Notebook retains the setup and associated result figures |
| 4.5 | `Nz = 40`, 30-m minimum depth, local `kz = 0.1 m2 s-1` | Day 12 minimum `12.72966` |
| 4.6 | `Nz = 100`, otherwise matching Trial 4.5 | Day 12 minimum `12.83041` |
| 4.8 | Strong local `kz = 1.0 m2 s-1` | Day 12 minimum `1.011435` |
| 4.9 | Depth-aware river spreading with local mixing | Day 12 minimum `23.42099` |
| 4.10 | Original river forcing with local `kz = 1e-3 m2 s-1` in the upper 30 m | Day 12 range `17.40515` to `40.14038` |
| 4.11 | Trial 4.10 plus ten-day ECCO salinity restoring | Day 12 range `18.90734` to `38.27020` |

Trials 4.7 and 4.12 are omitted because the saved notebooks do not contain a
completed result. Trial 4.0 and earlier superseded combined notebooks remain
available in Git history.

ECCO credentials are not stored in these notebooks. Set `ECCO_USERNAME` and
`ECCO_WEBDAV_PASSWORD` in the environment before running any trial that needs
to download ECCO data.
