# first-letters-clean-sheets

90 self-traced, mesh-level "clean" single-sheet patches of Herculaneum scroll surface (193.5 cm2 in all, 1.0-8.4 cm2 each, six scrolls, all scanned at 9.362 um), cut from traces grown with villa's own `vc_grow_seg_from_seed` on the public m7 surface prediction, each with a QA picture (flattened scan render, two CT cross-sections, frame-correct Reader v2 / d9v4C map at 40 um per pixel), plus the recipe for making more of them on a laptop, and a rough tile score ("on-sheet score") with its ranking of 120 public surfaces. The patches are meshes that were checked for one frame side, no folds and no stretched edges; the CT cross-sections show that most of them look like the flat public surfaces and unlike the swirled ones, but they do **not** prove that a patch lies on one papyrus layer. Reader v2 and d9v4C give speckle and blob fields on all 193.5 cm2 in both reading orders; we found no letters.

Status: contributed for the Vesuvius Challenge Progress Prize (October 2026). Everything here was made on one laptop with free compute. This is a local staging copy, not yet published.

## What is here

| path | what |
|---|---|
| `data/<scroll>/<patch>.tifxyz/` | the patch as a tifxyz (`x.tif`, `y.tif`, `z.tif`, `meta.json`), cropped to the clean rectangle. Coordinates are level-0 voxels (x, y, z) of the scan named in `meta.json`; `scale` 0.05 means nodes are about 20 voxels apart; `meta.json` has `voxel_size_um` (9.362), `volume_id`, `volume_shape_zyx`, the tracer `seed`, and a `clean_patch` block. |
| `index.csv` | one row per patch: scroll, volume id, voxel um, seed, area cm2, frame minority, fold share, frame-correct reading order, source run, and more (see below). |
| `qa/<scroll>/<patch>_qa.png` | the QA picture. `<patch>_both_orders.png` is the project's original ink-map picture in both reading orders. |
| `xsec/` | the CT cross-section measurements for every patch (10 sections each) and for six public control surfaces, plus the two gate checks. |
| `onsheet/` | the on-sheet tile score code (as run), thresholds, the 120-surface ranking table, group statistics and six example pictures. |
| `recipe/` | how the patches were made: Docker tracer on Windows, seed rule, mirror download, cutting a patch out of a trace. |
| `code/` | the audit, cross-section and picture code used here (`as_run/` holds the exact project scripts with their local paths). |
| `checks/` | recut reproduction, overlap between patches, integrity. |

Patches by scroll (from `index.csv`):

| scroll | volume id | patches | area cm2 | made in |
|---|---|---|---|---|
| PHerc0826 | 20250821151701 | 64 | 152.0 | 1 hand-seeded trace (8.4 cm2) + 30 + 33 grid-seeded |
| PHerc0813 | 20250821151723 | 15 | 24.6 | 15 |
| PHerc1545 | 20250821151648 | 4 | 7.7 | 4 |
| PHerc0358 | 20250821151737 | 5 | 6.3 | 5 |
| PHerc0125 | 20250821151825 | 1 | 1.7 | 1 |
| PHerc0257 | 20250821151750 | 1 | 1.2 | 1 |
| total | | 90 | 193.5 | |

A tracing run is still going on the same machine; later patches are not included in this snapshot.

`index.csv` columns beyond the ones above: node counts, `frame_side` (+1 normal points away from the scroll axis, -1 toward it), `rect_*` (checks recomputed on the cropped file), `trace_*` (the same checks on the whole trace the rectangle was cut from; the whole traces are often swirled, which is why only a rectangle is kept), `on_prediction_share`, the cross-section numbers, shares of the ink map above 0.5 / 0.8, and `flags` (listed under limits).

## How a patch was made and checked

1. **Trace.** Seeds are placed by a fixed rule on the public m7 surface prediction (one planar-band seed per 384-voxel cell, flattest first; `recipe/make_seeds.py`). Each seed is grown 75 generations with villa's `vc_grow_seg_from_seed` in the Docker image `ghcr.io/scrollprize/villa/volume-cartographer:main` (one container, four tracers side by side). The very first trace (PHerc0826, 8.4 cm2) was grown from a hand-picked seed on an official mesh. About 200 traces were run in all; 90 gave a kept patch.
2. **Mesh-level audit** (`code/tifxyz_audit.py`). A node is clean if it is finite, its normal differs from its neighbours' by at most 30 degrees, no edge is longer than 1.5 times the median edge, and it is on one frame side. The frame side is whether the normal points away from or toward the scroll axis (assumed to pass through the XY centre of the scan volume: a crude rule). The patch is the largest axis-aligned rectangle of clean nodes, with ground that earlier patches already cover (within 12 voxels) removed, and is kept if it is at least 1.0 cm2. `recipe/cut_patches.py` reproduces the stored rectangle exactly on all four stored traces we re-cut (`checks/`). No two patches of one scroll have any node within 12 voxels of each other (546,408 nodes tested).
3. **Reading.** Each patch was rendered (28 planes along the normal, 9.4 um per pixel) and read with Reader v2 and d9v4C in both reading orders. The frame-correct order (outward patch: F, inward patch: R) is the one that matches the scroll's winding at that patch under the project's rule. Maps are shown as min(Reader v2, d9v4C) at 40 um per pixel with a 1 cm bar. Model output is not ink and not letters.
4. **QA picture** (`qa/`): A, the flattened scan render (CT on the traced surface; fibres show as a crossed horizontal/vertical grid where the trace follows one sheet, and as swirled, streaked texture where it does not); B and C, two straightened CT cross-sections through the middle of the patch (depth -32..+32 voxels, about +-0.3 mm, along the normal, up = away from the scroll axis, the traced surface is the middle row, marked by red triangles); D, the frame-correct ink map.

### What the CT cross-sections show, and what they do not

Ten sections per patch (five along each direction) were measured and two are drawn. The measure is the tilt of the dominant layer orientation in the section (structure tensor, 0 degrees = layers parallel to the traced surface) and the drift this implies over the section length. It was fixed, with its gate checks and predictions, before any section was computed (private lab record X-218, 2026-10-19). Gates: the sampler reproduces the project renderer's stored stack to within one grey level at 100% of 1.8 million compared pixels; a synthetic shear of +-2 and +-4 degrees is recovered with median error 0.01 degrees (98% within 2 degrees, 5,960 windows); on public control surfaces the swirled ones (PHerc0826_v3_s0, PHerc0358_v3_s7, PHerc0257_v2_s0, chosen from our earlier visual observation) show a median |tilt| of 4.70 degrees against 1.01 degrees for the flat ones (PHerc0826_auto_grown_...841364, ...849408, ...900360), a ratio of 4.6.

Result on the 90 patches (900 sections): median |tilt| 0.94 degrees; per-patch medians 0.5-1.6 degrees at the 10th-90th percentile, up to 3.6 degrees. That is the same as the flat controls and about five times below the swirled ones. Most patch medians lie in the flat-control range; the highest one (3.6 degrees) is within the range of the swirled controls' per-surface medians (3.2-7.7 degrees). But:

- Over a section that is 1,000-4,000 voxels long, 1 degree is already 17-70 voxels of drift, and the layers in these sections are only about 10-12 voxels apart (see the pictures: a stack of wavy laminations fills the whole +-32 voxels). We had set 20 voxels as the line for "about one layer spacing"; 57% of patch sections exceed it, and so do 77% of the flat control sections and 90% of the swirled ones. Our advance expectation that at least 80% of patches would stay under it was not met. So this measure rules out gross swirl, and it cannot confirm that a patch stays on one layer.
- By eye, the sections show laminations across the whole depth, so a section cannot tell "this sheet" from its neighbours; some patches show a clear horizontal band at the traced surface (for example `PHerc0826_t269_s051`), others, such as `PHerc0813_t269_s000`, show a regular fibre grid on one half and swirled texture on the other.
- Six patches have a median |tilt| of 2 degrees or more (twice the flat control median; a post-hoc cut, flagged `xsec_tilt_ge_2deg_post_hoc`): PHerc0125_t268_s029 (3.6, the only PHerc0125 patch), PHerc0826_t268_s012, s019, s033, PHerc0826B_t269_s016, PHerc0813_t269_s000.
- The "informative section" rule (coherence above a pixel-shuffle noise null) was true for all 900 sections, so it does not screen anything.
- The cross-sections use the same trilinear sampler as the renders (agreement with the stored stack is about software, not about the sheet). Cached CT was used; 0.35% of section columns fell in chunks we did not have and were marked invalid.

Also, post-hoc: the share of each patch's nodes that lie on the public m7 prediction it was grown on has a median 0.52 (10th percentile 0.41, minimum 0.15; shifting the same nodes 15-40 voxels along the normal gives about 0.25). Five patches have less than 0.30 or were named in an earlier note (`off_public_prediction`): PHerc0826_t268_s004, s014, s033, PHerc0826_t269_s052, s058.

## Ink maps

On all 90 patches, in the frame-correct order, the share of pixels with min(Reader v2, d9v4C) above 0.5 has median 0.042 (maximum 0.49) and the share above 0.8 is 0 in every patch (the other order is similar). Looks at the pictures by eye (three looks, the first not blind, the other two blind to the scores) found blob and speckle fields, no rows of letter-shaped forms. Two features that looked interesting were followed up and were not letter-like: a short band in `PHerc0813_t269_s012` (reversed order only) was not reproduced as a band by four other readers; a column of blobs in `PHerc0826_t269_s101` (frame-correct order) is reproduced as blobs by two of four other readers (weakly by the other two) and does not persist when the surface is shifted 17 voxels off the sheet. Five patches have a few non-finite reader pixels (flag `ink_map_nonfinite_px`, 1-2,086 pixels of about 1.4 million). We have no known-ink control on these scrolls, and both readers were not trained for this regime, so "no letters found" is a statement about these model maps over these patches, not about whether writing is absent.

## Where tracing worked and where it did not

The patch yield depends strongly on where the tracer is started. Blocks chosen by the densest public surface prediction (PHerc0125, PHerc0257) gave 2.9 cm2 from 32 traces; a paired control block on PHerc0813 chosen the same way gave 0 cm2 from 16 traces. Blocks chosen around the mesh that the on-sheet score ranks best for that scroll gave 24.6 cm2 from 24 traces on PHerc0813 (a number still growing when this snapshot was taken), 6.3 cm2 from 16 on PHerc0358, 7.7 cm2 from 16 on PHerc1545, and PHerc0826 gave 75.8, 40.7 and 27.1 cm2 (the last a second box) from three runs of 32 traces. The cause of the difference is not determined: thin versus fused mask thickness does not explain it. Each of these comparisons is a single experiment, not a replicated one.

## The on-sheet score and the 120-surface ranking

`onsheet/` contains the code (as run; it needs 28-plane renders of surfaces, so it is documentation more than a tool), the thresholds (taken only from a known-writing surface, frozen before any public surface was scored) and `ranking_120_surfaces.csv`: 120 public third-party surfaces (46 `auto_grown` meshes from flummoxjr/gp13-ink-detectability, 50 from the Bullo27 first-letters-survey, 20 from kadenpool/scroll-lineup, 4 from pscamillo/vesuvius-eligible-meshes) ranked by on-sheet area. The score marks a 256-pixel tile as on-sheet if its fibre texture is strongly oriented, its depth contrast is within the range seen on known writing, and its texture is above a floor. It is a rough screen: held-out known writing scores 0.74 / 0.74 / 0.96 on-sheet, shuffled nulls 0.000 and 0.005, but on the six surfaces that we had looked at by eye (two flat, four swirled) it separates them with AUC 0.875 only, and one clearly swirled surface (PHerc0826_v3_s0, share 0.35) scores above a flat one (PHerc0800, 0.14).

What it gives: median on-sheet share 0.19 over the 120 surfaces; 15 of 120 reach 0.5; 199.5 of 1,007.5 cm2 pass. The `auto_grown` meshes score highest (median share 0.36, 14 of 15 high scorers), the survey meshes lowest (median 0.08, 1 of 50 at 0.5 or more).

The finding we can stand behind is therefore narrower than "public auto-grown surfaces cross layers". On the few contact sheets we inspected (`onsheet/examples/` shows the six surfaces of observation O-347 in our records), the scan texture of most surfaces was swirled and layer-crossing rather than a regular fibre grid, and the regular-grid ones were a minority; this was a visual count on a handful of sheets, not a count over all 120, and the swirled examples were mostly survey-type surfaces (v2/v3 segments), not the `auto_grown` meshes. The tile score does not detect swirl reliably, and the cross-section measure above (tested on three of the swirled and three flat surfaces) separates them better. Surfaces made by others for other purposes (ink hunting, tests) are not criticised here; we use them as controls and a baseline only, and their meshes are not redistributed.

## Honest limits

- Small patches (1.0-8.4 cm2, median 1.8) and mesh-level checks: a clean mesh is not shown to lie on one papyrus layer (see the cross-section section); swirled ground can pass the mesh test, and 6 of the 24 official PHerc0826 meshes also pass a 4 cm2 rectangle.
- Agreement of a trace with the public prediction is self-consistency (same tool, same prediction), not independent confirmation of the sheet. The extent of a trace beyond its first 75 generations has lower prediction support.
- Traces are not bit-exact (two runs of one seed differ by about 2.4 voxels median); the side convention of the output depends on the seed, so the frame has to be measured on every trace.
- The frame rule uses the XY centre of the scan volume as the scroll axis.
- Ink maps are model output at the edge of what the readers were built for; no known-ink control on these scrolls; no letters were found.
- Several patches are flagged in `index.csv` (5 `off_public_prediction`, 5 `ink_map_nonfinite_px`, 3 `trace_touches_mirror_box_wall`, 6 `xsec_tilt_ge_2deg_post_hoc`); they are included, not dropped.
- Cropping changes border normals: 21 patches have a few fold nodes within the outer two node rings of the cropped file (interior fold share is 0 everywhere); `rect_fold_share_interior` and `checks/` say so.
- Areas: `area_cm2` is nodes times (median edge x voxel size)^2 as in the original audit; the triangulated surface area of the crop is 90-100% of it.
- Made with villa's own tracer and the public prediction: this is not a new tracing method.

## How to reproduce

- Audit a patch: `python code/tifxyz_audit.py data/PHerc0813/PHerc0813_t269_s000.tifxyz --shape 16993 7947 7947` (needs numpy, tifffile).
- Cross-sections: `python code/xsec_patch.py <patch.tifxyz> --ct <dir of level-0 raw 128^3 chunks>` (numpy, scipy, tifffile; the public level-0 zarr stores these chunks uncompressed; fetch only the chunks you need, a patch needs a few hundred MB). Re-running on `PHerc0813_t269_s014` reproduced our stored tilts to the printed digit.
- Make more patches: `recipe/README.md` (download a mirror block, Docker volume, trace rounds, cut). Needs about 6 GB of download and a few hours per 16 traces; Windows with Docker Desktop was used. The seed rule (`make_seeds.py`) reproduced all 9 stored seeds we tried; `fetch_mirror.py` and `trace_round.py` are cleaned copies of code that ran and were not re-run in this form.
- The QA pictures need the stored 28-plane stacks and reader maps, which are not included (about 6 GB; the readers' weights are public); `code/as_run/` shows how they were built.

## Credits

Vesuvius Challenge: the scans, the open-data bucket, the public m7 surface prediction and its normal grids, and villa's volume-cartographer (`vc_grow_seg_from_seed`, `ghcr.io/scrollprize/villa/volume-cartographer:main`, github.com/ScrollPrize/villa). Reader v2: DomRusso2/reader-v2. d9v4C: TAUIL's shape-loss fine-tune of Reader v2. Survey and surface lists we used as controls and baselines: Bullo27 (first-letters-survey), pscamillo (vesuvius-eligible-meshes and his published ink maps), flummoxjr (gp13-ink-detectability), kadenpool (scroll-lineup). We did not write any of these tools or models; the patches, checks, cross-section measure and recipe are ours.

## Licence

Code (`code/`, `recipe/`, `onsheet/code_as_run/`, `checks/` scripts): MIT, see `LICENSE`. Derived data (the tifxyz patches, pictures, tables): CC BY-NC 4.0, like the Vesuvius Challenge scans they derive from; see `LICENSE-DATA.md`.
