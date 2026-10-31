# first-letters-clean-sheets

173 self-traced, mesh-level "clean" single-sheet patches of Herculaneum scroll surface (355.0 cm2 in all, 1.0-8.4 cm2 each, on 13 First Letters-eligible scrolls, scanned at 9.362 um or 8.64 um), cut from traces grown with villa's own `vc_grow_seg_from_seed` on the public m7 surface prediction, each with a QA picture (flattened scan render, two CT cross-sections, frame-correct Reader v2 / d9v4C map at 40 um per pixel), plus the recipe for making more of them on a laptop, and a rough tile score ("on-sheet score") with its ranking of 120 public surfaces. The patches are meshes that were checked for one frame side, no folds and no stretched edges; the CT cross-sections show that most of them look like the flat public surfaces and unlike the swirled ones, but they do **not** prove that a patch lies on one papyrus layer. Reader v2 and d9v4C give speckle and blob fields on all 355.0 cm2 in both reading orders; we found no letters.

Status: contributed for the Vesuvius Challenge Progress Prize (October 2026). Everything here was made on one laptop with free compute. Published 2026-10-19 with 90 patches; updated 2026-10-21 with 83 more from the finished tracing runs (seven of the 13 scrolls are new in this update).

## What is here

| path | what |
|---|---|
| `data/<scroll>/<patch>.tifxyz/` | the patch as a tifxyz (`x.tif`, `y.tif`, `z.tif`, `meta.json`), cropped to the clean rectangle. Coordinates are level-0 voxels (x, y, z) of the scan named in `meta.json`; `scale` 0.05 means nodes are about 20 voxels apart; `meta.json` has `voxel_size_um` (9.362 or 8.64), `volume_id`, `volume_shape_zyx`, the tracer `seed`, and a `clean_patch` block. |
| `index.csv` | one row per patch: scroll, volume id, voxel um, seed, area cm2, frame minority, fold share, frame-correct reading order, source run, and more (see below). The five `lh_*` columns are LynnHenry's independent sheet-contrast check and reading order for `vc_render_tifxyz` default renders (without `--flip-normals`) for each patch (from [LynnHenry/sheet-contrast-check](https://github.com/LynnHenry/sheet-contrast-check) `data/clean_sheets_pass_fail.csv`, commit 53486e7): `lh_sheet_check` is pass for the 79 patches whose sheet contrast reaches the weakest curated control of the same scan protocol. |
| `qa/<scroll>/<patch>_qa.png` | the QA picture. `<patch>_both_orders.png` is the project's original ink-map picture in both reading orders. |
| `xsec/` | the CT cross-section measurements for every patch (10 sections each) and for six public control surfaces, plus the two gate checks. |
| `onsheet/` | the on-sheet tile score code (as run), thresholds, the 120-surface ranking table, group statistics and six example pictures. |
| `recipe/` | how the patches were made: Docker tracer on Windows, seed rule, mirror download, cutting a patch out of a trace. |
| `code/` | the audit, cross-section and picture code used here (`as_run/` holds the exact project scripts with their local paths). |
| `checks/` | recut reproduction, overlap between patches, integrity. |

Patches by scroll (from `index.csv`):

| scroll | volume id | scan um | patches | area cm2 | made in |
|---|---|---|---|---|---|
| PHerc0826 | 20250821151701 | 9.362 | 80 | 185.9 | 1 hand-seeded trace (8.4 cm2) + 30 + 33 + 16 grid-seeded (the last 16 in a continuation of the second box) |
| PHerc0813 | 20250821151723 | 9.362 | 21 | 37.9 | 17 + 4 (continuation) |
| PHerc0175A | 20250521115057 | 8.64 | 15 | 32.7 | 15 (new in the update) |
| PHerc0211 | 20250821151803 | 9.362 | 11 | 23.6 | 11 (new) |
| PHerc1545 | 20250821151648 | 9.362 | 13 | 23.5 | 4 + 9 (continuation) |
| PHerc1203 | 20250820131727 | 9.362 | 9 | 14.3 | 9 (new) |
| PHerc0306B | 20250521133212 | 8.64 | 7 | 11.6 | 7 (new) |
| PHerc0191 | 20250821151635 | 9.362 | 5 | 7.5 | 5 (new) |
| PHerc0483B | 20251124083638 | 8.64 | 3 | 6.6 | 3 (new) |
| PHerc0358 | 20250821151737 | 9.362 | 5 | 6.3 | 5 |
| PHerc0846B | 20250804142305 | 9.362 | 2 | 2.3 | 2 (new) |
| PHerc0125 | 20250821151825 | 9.362 | 1 | 1.7 | 1 |
| PHerc0257 | 20250821151750 | 9.362 | 1 | 1.2 | 1 |
| total | | | 173 | 355.0 | |

Every volume above is on the First Letters eligible list as published on scrollprize.org/prizes (checked 2026-10-21, 21 volumes). Two PHerc0813 patches (`PHerc0813_t269_s024`, `s040`, 3.3 cm2) come from the first run but finished after the first snapshot. `PHerc0826B_t275_s064` (2.8 cm2) was rendered and read later than the rest: its first render failed on a lock-file error during the run and was redone afterwards with the same code. The tracing runs have ended; 16 traces on an eighth new scroll, PHerc0490B, gave no patch.

`index.csv` columns beyond the ones above: node counts, `frame_side` (+1 normal points away from the scroll axis, -1 toward it), `rect_*` (checks recomputed on the cropped file), `trace_*` (the same checks on the whole trace the rectangle was cut from; the whole traces are often swirled, which is why only a rectangle is kept), `on_prediction_share`, the cross-section numbers, shares of the ink map above 0.5 / 0.8, and `flags` (listed under limits).

## How a patch was made and checked

1. **Trace.** Seeds are placed by a fixed rule on the public m7 surface prediction (one planar-band seed per 384-voxel cell, flattest first; `recipe/make_seeds.py`). Each seed is grown 75 generations with villa's `vc_grow_seg_from_seed` in the Docker image `ghcr.io/scrollprize/villa/volume-cartographer:main` (one container, four tracers side by side). The very first trace (PHerc0826, 8.4 cm2) was grown from a hand-picked seed on an official mesh. About 400 traces were run in all; 173 gave a kept patch, all of which are here.
2. **Mesh-level audit** (`code/tifxyz_audit.py`). A node is clean if it is finite, its normal differs from its neighbours' by at most 30 degrees, no edge is longer than 1.5 times the median edge, and it is on one frame side. The frame side is whether the normal points away from or toward the scroll axis (assumed to pass through the XY centre of the scan volume: a crude rule). The patch is the largest axis-aligned rectangle of clean nodes, with ground that earlier patches already cover (within 12 voxels) removed, and is kept if it is at least 1.0 cm2. `recipe/cut_patches.py` reproduces the stored rectangle exactly on all four stored traces we re-cut (`checks/`). No two patches of one scroll have any node within 12 voxels of each other (1,017,802 nodes tested).
3. **Reading.** Each patch was rendered (28 planes along the normal, 9.4 um per pixel) and read with Reader v2 and d9v4C in both reading orders. The frame-correct order (outward patch: F, inward patch: R) is the one that matches the scroll's winding at that patch under the project's rule. **F/R and `frame_correct_order` refer to the layer order of our own renders (the project renderer used for the QA pictures), not to `vc_render_tifxyz` output.** LynnHenry measured that `vc_render_tifxyz` stacks the layers in the opposite order, and that for `ink_9um`-family models on `vc_render_tifxyz` renders the right order is the reverse of our field (outward patch: reverse); see [LynnHenry/sheet-contrast-check](https://github.com/LynnHenry/sheet-contrast-check) for the calibration on labelled PHerc0841. We checked this ourselves on three patches (PHerc0175A_t275_s000, PHerc0813_t269_s002, PHerc0826_t268_s000) with villa's `vc_render_tifxyz` (VC3D build d57c28e, 28 slices, CT streamed from the open-data bucket). With `--flip-normals`, which Tutorial 5 uses, the layers come out in our order (plane-by-plane r 0.975-0.987 with our stacks), so our F/R apply as they are. With the tool's default, the layers are reversed, so flip our field. That is what the `lh_order_for_vc_render_tifxyz` column gives, and it matched on all three. Reader V3 run through villa's own flat inference on these renders reproduces our maps (r 0.88-0.95 in the matching order). The same repository reports a sheet-contrast check on all 173 patches (79 with a sheet signal in the range of curated controls), which is a useful complement to our cross-section tilt check. Maps are shown as min(Reader v2, d9v4C) at 40 um per pixel with a 1 cm bar. Model output is not ink and not letters.
4. **QA picture** (`qa/`): A, the flattened scan render (CT on the traced surface; fibres show as a crossed horizontal/vertical grid where the trace follows one sheet, and as swirled, streaked texture where it does not); B and C, two straightened CT cross-sections through the middle of the patch (depth -32..+32 voxels, about +-0.3 mm, along the normal, up = away from the scroll axis, the traced surface is the middle row, marked by red triangles); D, the frame-correct ink map.

### What the CT cross-sections show, and what they do not

Ten sections per patch (five along each direction) were measured and two are drawn. The measure is the tilt of the dominant layer orientation in the section (structure tensor, 0 degrees = layers parallel to the traced surface) and the drift this implies over the section length. It was fixed, with its gate checks and predictions, before any section was computed (private lab record X-218, 2026-10-19). Gates: the sampler reproduces the project renderer's stored stack to within one grey level at 100% of 1.8 million compared pixels; a synthetic shear of +-2 and +-4 degrees is recovered with median error 0.01 degrees (98% within 2 degrees, 5,960 windows); on public control surfaces the swirled ones (PHerc0826_v3_s0, PHerc0358_v3_s7, PHerc0257_v2_s0, chosen from our earlier visual observation) show a median |tilt| of 4.70 degrees against 1.01 degrees for the flat ones (PHerc0826_auto_grown_...841364, ...849408, ...900360), a ratio of 4.6.

Result on the 173 patches (1,730 sections): median |tilt| 0.97 degrees (0.94 on the first 90, 1.00 on the 83 added); per-patch medians 0.5-1.8 degrees at the 10th-90th percentile, up to 4.0 degrees. That is the same as the flat controls and about five times below the swirled ones. Most patch medians (131 of 173) lie at or below the highest flat-control median (1.3 degrees); the two highest (4.0 and 3.6 degrees) are within the range of the swirled controls' per-surface medians (3.2-7.7 degrees). The cross-section measure, its parameters and the controls were not changed for the added patches. But:

- Over a section that is 1,000-4,000 voxels long, 1 degree is already 17-70 voxels of drift, and the layers in these sections are only about 10-12 voxels apart (see the pictures: a stack of wavy laminations fills the whole +-32 voxels). We had set 20 voxels as the line for "about one layer spacing"; 58% of patch sections exceed it, and so do 77% of the flat control sections and 90% of the swirled ones. Our advance expectation that at least 80% of patches would stay under it was not met. So this measure rules out gross swirl, and it cannot confirm that a patch stays on one layer.
- By eye, the sections show laminations across the whole depth, so a section cannot tell "this sheet" from its neighbours; some patches show a clear horizontal band at the traced surface (for example `PHerc0826_t269_s051`), others, such as `PHerc0813_t269_s000`, show a regular fibre grid on one half and swirled texture on the other.
- 13 patches have a median |tilt| of 2 degrees or more (twice the flat control median; a post-hoc cut, flagged `xsec_tilt_ge_2deg_post_hoc`): PHerc0846B_t275_s008 (4.0, one of the two PHerc0846B patches; the other, s009, is 2.0), PHerc0125_t268_s029 (3.6, the only PHerc0125 patch), PHerc1203_t275_s004 (3.0), PHerc0826_t268_s012, s019, s033, PHerc0826B_t269_s016, PHerc0813_t269_s000, PHerc0191_t275_s001, s007, PHerc0483B_t275_s005, PHerc1545_t275_s054 (2.0-2.3 each).
- The "informative section" rule (coherence above a pixel-shuffle noise null) was true for all 1,720 sections, so it does not screen anything.
- The cross-sections use the same trilinear sampler as the renders (agreement with the stored stack is about software, not about the sheet). Cached CT was used; 0.36% of section columns fell in chunks we did not have and were marked invalid.

Also: the share of each patch's nodes that lie on the public m7 prediction it was grown on has a median 0.67 (10th percentile 0.46, minimum 0.15; shifting the same nodes 15-40 voxels along the normal gives about 0.25). For the 80 patches of the later run this check was planned in advance (all at 0.44 or more); for the others it is post-hoc. Five patches have less than 0.30 or were named in an earlier note (`off_public_prediction`): PHerc0826_t268_s004, s014, s033, PHerc0826_t269_s052, s058.

## Ink maps

On all 173 patches, in the frame-correct order, the share of pixels with min(Reader v2, d9v4C) above 0.5 has median 0.051 (maximum 0.49) and the share above 0.8 is 0 in every patch (the other order is similar). Looks at the pictures of the first 90 by eye (three looks, the first not blind, the other two blind to the scores) found blob and speckle fields, no rows of letter-shaped forms. The 83 added patches were looked at once by eye, not blind (the later run's per-scroll contact sheets and the two late PHerc0813 pictures and the late `PHerc0826B_t275_s064` picture): the same blob and speckle fields on all seven new scrolls, no rows of letter-shaped forms. Two features that looked interesting were followed up and were not letter-like: a short band in `PHerc0813_t269_s012` (reversed order only) was not reproduced as a band by four other readers; a column of blobs in `PHerc0826_t269_s101` (frame-correct order) is reproduced as blobs by two of four other readers (weakly by the other two) and does not persist when the surface is shifted 17 voxels off the sheet. Seven patches have non-finite reader pixels (flag `ink_map_nonfinite_px`): six have 1-2,086 pixels; `PHerc0826B_t275_s050` (1.05 cm2) has 49,302, a jagged strip along one edge of its map (cause not determined; the patch is not at the edge of the scan volume). We have no known-ink control on these scrolls, and both readers were not trained for this regime, so "no letters found" is a statement about these model maps over these patches, not about whether writing is absent.

## Where tracing worked and where it did not

The patch yield depends strongly on where the tracer is started. Blocks chosen by the densest public surface prediction (PHerc0125, PHerc0257) gave 2.9 cm2 from 32 traces; a paired control block on PHerc0813 chosen the same way gave 0 cm2 from 16 traces. Blocks chosen around the mesh that the on-sheet score ranks best for that scroll gave 24.6 cm2 from the first 24 traces on PHerc0813 (27.9 cm2 from all 32 of that run), 6.3 cm2 from 16 on PHerc0358, 7.7 cm2 from 16 on PHerc1545, and PHerc0826 gave 75.8, 40.7 and 27.1 cm2 (the last a second box) from three runs of 32 traces. The cause of the difference is not determined: thin versus fused mask thickness does not explain it. Each of these comparisons is a single experiment, not a replicated one, and across the five anchored blocks the on-sheet score of the anchor predicts the yield only weakly (Spearman 0.30).

The later run (2026-10-19 to 10-21) applied the same block rule, unchanged, to eight more eligible scrolls (16 traces each) and continued three earlier blocks with later, less flat seeds. New blocks: PHerc0175A 32.7 cm2, PHerc0211 23.6, PHerc1203 14.3, PHerc0306B 11.6, PHerc0191 7.5, PHerc0483B 6.6, PHerc0846B 2.3, PHerc0490B 0. PHerc0175A is the only new block where the whole traces were mostly on one frame side (median off-side share 0.07, against 0.12-0.41 on the others). Continued blocks: PHerc0813 +10.0 cm2 from 16 more traces, PHerc0826 second box +33.9 from 32 (31.1 included), PHerc1545 +15.8 from 24; many later seeds there sat on sheet already traced and were skipped. These too are single runs.

## The on-sheet score and the 120-surface ranking

`onsheet/` contains the code (as run; it needs 28-plane renders of surfaces, so it is documentation more than a tool), the thresholds (taken only from a known-writing surface, frozen before any public surface was scored) and `ranking_120_surfaces.csv`: 120 public third-party surfaces (46 `auto_grown` meshes from flummoxjr/gp13-ink-detectability, 50 from the Bullo27 first-letters-survey, 20 from kadenpool/scroll-lineup, 4 from pscamillo/vesuvius-eligible-meshes) ranked by on-sheet area. The score marks a 256-pixel tile as on-sheet if its fibre texture is strongly oriented, its depth contrast is within the range seen on known writing, and its texture is above a floor. It is a rough screen: held-out known writing scores 0.74 / 0.74 / 0.96 on-sheet, shuffled nulls 0.000 and 0.005, but on the six surfaces that we had looked at by eye (two flat, four swirled) it separates them with AUC 0.875 only, and one clearly swirled surface (PHerc0826_v3_s0, share 0.35) scores above a flat one (PHerc0800, 0.14).

What it gives: median on-sheet share 0.19 over the 120 surfaces; 15 of 120 reach 0.5; 199.5 of 1,007.5 cm2 pass. The `auto_grown` meshes score highest (median share 0.36, 14 of 15 high scorers), the survey meshes lowest (median 0.08, 1 of 50 at 0.5 or more).

The finding we can stand behind is therefore narrower than "public auto-grown surfaces cross layers". On the few contact sheets we inspected (`onsheet/examples/` shows the six surfaces of observation O-347 in our records), the scan texture of most surfaces was swirled and layer-crossing rather than a regular fibre grid, and the regular-grid ones were a minority; this was a visual count on a handful of sheets, not a count over all 120, and the swirled examples were mostly survey-type surfaces (v2/v3 segments), not the `auto_grown` meshes. The tile score does not detect swirl reliably, and the cross-section measure above (tested on three of the swirled and three flat surfaces) separates them better. Surfaces made by others for other purposes (ink hunting, tests) are not criticised here; we use them as controls and a baseline only, and their meshes are not redistributed.

## Honest limits

- Small patches (1.0-8.4 cm2, median 1.7) and mesh-level checks: a clean mesh is not shown to lie on one papyrus layer (see the cross-section section); swirled ground can pass the mesh test, and 6 of the 24 official PHerc0826 meshes also pass a 4 cm2 rectangle.
- Agreement of a trace with the public prediction is self-consistency (same tool, same prediction), not independent confirmation of the sheet. The extent of a trace beyond its first 75 generations has lower prediction support.
- Traces are not bit-exact (two runs of one seed differ by about 2.4 voxels median); the side convention of the output depends on the seed, so the frame has to be measured on every trace.
- The frame rule uses the XY centre of the scan volume as the scroll axis.
- Ink maps are model output at the edge of what the readers were built for; no known-ink control on these scrolls; no letters were found.
- Several patches are flagged in `index.csv` (5 `off_public_prediction`, 7 `ink_map_nonfinite_px`, 3 `trace_touches_mirror_box_wall`, 13 `xsec_tilt_ge_2deg_post_hoc`); they are included, not dropped.
- Cropping changes border normals: 51 patches have a few fold nodes within the outer two node rings of the cropped file (interior fold share is 0 everywhere); `rect_fold_share_interior` and `checks/` say so.
- Areas: `area_cm2` is nodes times (median edge x voxel size)^2 as in the original audit; the triangulated surface area of the crop is 87-105% of it.
- Made with villa's own tracer and the public prediction: this is not a new tracing method.

## How to reproduce

- Audit a patch: `python code/tifxyz_audit.py data/PHerc0813/PHerc0813_t269_s000.tifxyz --shape 16993 7947 7947` (needs numpy, tifffile).
- Cross-sections: `python code/xsec_patch.py <patch.tifxyz> --ct <dir of level-0 raw 128^3 chunks>` (numpy, scipy, tifffile; the public level-0 zarr stores these chunks uncompressed; fetch only the chunks you need, a patch needs a few hundred MB). Re-running on `PHerc0813_t269_s014` reproduced our stored tilts to the printed digit.
- Make more patches: `recipe/README.md` (download a mirror block, Docker volume, trace rounds, cut). Needs about 6 GB of download and a few hours per 16 traces; Windows with Docker Desktop was used. The seed rule (`make_seeds.py`) reproduced all 9 stored seeds we tried; `fetch_mirror.py` and `trace_round.py` are cleaned copies of code that ran and were not re-run in this form.
- The QA pictures need the stored 28-plane stacks and reader maps, which are not included (about 6 GB; the readers' weights are public); `code/as_run/` shows how they were built.

## Authors

Jack Gilman and Andrew Gilman.

## Credits

Vesuvius Challenge: the scans, the open-data bucket, the public m7 surface prediction and its normal grids, and villa's volume-cartographer (`vc_grow_seg_from_seed`, `ghcr.io/scrollprize/villa/volume-cartographer:main`, github.com/ScrollPrize/villa). Reader v2: DomRusso2/reader-v2. d9v4C: TAUIL's shape-loss fine-tune of Reader v2. Survey and surface lists we used as controls and baselines: Bullo27 (first-letters-survey), pscamillo (vesuvius-eligible-meshes and his published ink maps), flummoxjr (gp13-ink-detectability), kadenpool (scroll-lineup). LynnHenry ([sheet-contrast-check](https://github.com/LynnHenry/sheet-contrast-check)): the per-patch sheet-contrast results and `vc_render_tifxyz` reading order in the `lh_*` columns of `index.csv` (CC BY-NC 4.0), and the finding that our F/R convention is reversed for `vc_render_tifxyz` renders. We did not write any of these tools or models; the patches, checks, cross-section measure and recipe are ours.

## Licence

Code (`code/`, `recipe/`, `onsheet/code_as_run/`, `checks/` scripts): MIT, see `LICENSE`. Derived data (the tifxyz patches, pictures, tables): CC BY-NC 4.0, like the Vesuvius Challenge scans they derive from; see `LICENSE-DATA.md`.
