# CT cross-section measurements (X-218)

`patch_sections.csv`: 10 sections per patch (columns: canvas row/column of the line, valid columns, arc length, tilt in degrees, coherence, noise threshold, informative flag (always 1), implied drift in voxels,
uncached columns). `control_sections.csv`: the same for 3 flat and 3 swirled public surfaces. `summary.json`: all aggregate numbers and per-patch medians. `g1.json`: synthetic-shear gate. `identity.json`: the
sampler-versus-stored-render check. Code: `code/xsec.py`, `code/xsec_patch.py`. Definition of tilt and drift: README, section "What the CT cross-sections show". Positive depth = away from the scroll axis.
