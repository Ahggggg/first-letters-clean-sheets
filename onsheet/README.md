# On-sheet score and ranking of 120 public surfaces

`ranking_120_surfaces.csv`: 120 public third-party surfaces ranked by on-sheet area (columns: rank, surface, scroll, source group, on-sheet area cm2, valid area cm2, on-sheet share, tile counts, rank by share, median
depth-contrast and texture). `ranking_stats.json`: group statistics. `thresholds259.json` (+ sha256): thresholds frozen from a known-writing surface (w00) before any of the 120 was scored.
`code_as_run/`: the scripts as run (project paths inside; they need 28-plane renders of each surface; `onsheet259.py` is the final rule, `onsheet257*.py` the tile code it imports, `pic245.py` a helper).

Rule (per 256 x 256 canvas-pixel tile, at least half valid): A = fourth-harmonic orientation coherence of the surface render (crossed fibres score high), C = range of the tile-mean intensity profile over the 28 planes
in percent of its mean (must lie between the known-writing 10th and 99th percentile), Tx = gradient texture of the surface render above the known-writing 10th percentile. On-sheet iff all three pass.

Validation we have (all in the project record X-202): held-out known writing 0.74 / 0.74 / 0.96 of tiles on-sheet; shuffled-pixel null 0.000, spatial-shuffle null 0.005; on six surfaces we had looked at by eye
(two flat, four swirled) AUC 0.875, and the unseen pair at 0.5; one swirled surface (PHerc0826_v3_s0, 0.35) scores above one flat surface (PHerc0800, 0.14). So: a rough screen, not a swirl detector.

`examples/`: ex1 (a regular fibre-grid surface, 95% on-sheet), ex2 (rank 3, scores high but is visibly swirled: a false positive), ex3-ex4 (flat references), ex5-ex6 (swirled). Each picture is scan plane 14 | outward mask |
writing-order map (min of the two readers, F where plane 0 faces outward, R where inward, grey where undetermined), 40 um per pixel. `onsheet_map_all_120.png` shows the on-sheet tiles of all 120.
By source group the median on-sheet share is: flummoxjr auto_grown 0.36 (n 46), kadenpool 0.17 (20), pscamillo 0.22 (4), Bullo27 survey 0.08 (50). The surfaces were made by others for other purposes; they are
baselines here, not targets of criticism, and their meshes are not redistributed.
