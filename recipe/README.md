# Recipe: tracing clean sheets with villa's tracer on one Windows laptop

What we did, in order, and what bit us. Hardware: 4-core laptop CPU (5 CPUs given to Docker Desktop, 14 GB), no GPU used for tracing. Time: about 15-20 minutes per tracer with four side by side
(75 generations), so one round of four traces is about 20 minutes, a 16-trace block about 1.5 hours.

1. **Pick a block (a box of a few thousand voxels).** The first lesson: anchor it on a mesh that already looks like a single sheet. Our blocks chosen by the densest public surface prediction gave almost nothing
   (PHerc0125 and PHerc0257: 2.9 cm2 from 32 traces; a PHerc0813 control chosen that way: 0 cm2 from 16), while blocks centred on the best-ranked mesh of the on-sheet score gave 24.6 cm2 from 24 traces
   (PHerc0813), 6.3 (PHerc0358), 7.7 (PHerc1545), 27.1 (a second PHerc0826 box). Rule used: take the anchor mesh's bounding box, widen it by 700 voxels (the seed wall margin) on every side, make each side at
   least 3,456 and at most 5,760 voxels, round outward to the 192-voxel chunk grid (`as_run/` in the project's `sites269.py`, not copied here; the idea is all that matters). One experiment, not replicated.
2. **Download the block** (`fetch_mirror.py SCROLL VOLUME_ID --box ... --out DIR`): level-0 chunks of the public m7 surface prediction inside the box (0.4-1.5 GB) and the normal-grid slices (3-5 GB, whole slices).
   Public bucket path: `<scroll>/representations/predictions/surfaces/<volume>-surface-20260413222639-surface-m7-L0-th0.2` (`.zarr` and `.normal-grids`).
3. **Put the mirror into a Docker named volume, not a bind mount** (`trace_round.py volume`). From a Windows bind mount the tracer sat in state D at about 0.3 of a core (25 generations in 28 minutes,
   aborted); from a named volume it is about 100 times faster (30 generations in 62 s on the small box). Copy of 3-6 GB: about 7 minutes. Docker volumes live on the system drive: keep 10+ GB free.
4. **The voxel size must be in `meta.json` at the volume root** (`{"type": "vol", "uuid", "width", "height", "slices", "voxelsize": 9.362}`). The public `metadata.json` has no `scan` key and makes the tracer
   throw; and without a voxel size the tracer computes area 0, applies `min_area_cm` (0.3) and silently DELETES the finished output. `trace_round.py volume` writes it. 9.362 um for PHerc0125, 0191, 0211, 0257,
   0358, 0813, 0826, 0846, 1203, 1545; the 8.64 um scrolls need 8.64.
5. **Seeds** (`make_seeds.py MIRROR --box ... --out seed-candidates.json`): one planar-band seed per 384-voxel cell, flattest first (rule in the file). Skip a seed within 16 voxels of an earlier whole trace
   (same sheet already grown) and keep the seeds of a round at least 1,000 voxels apart. Reproduced all 9 seeds of the PHerc0813 mirror that we compared (cell by cell).
6. **Trace** (`trace_round.py round --name ... --voxel-um 9.362 --out ... --seeds "x y z" ...`): four tracers in one detached container, parameters of the official PHerc0826 recipe: mode seed, generations 75,
   thread_limit 1, use_cuda false, min_area_cm 0.3, cache_size 3e9; entry point overridden so no AWS credentials or network are needed; container removed afterwards. Under 0.5 GB memory per tracer in our runs.
   (Run tracers detached and hidden if you want your laptop back: a Windows console window per container is annoying.)
7. **Audit and cut** (`cut_patches.py TRACE/seg --shape Z Y X --exclude earlier patches ... --out ...`): measure the frame side of every trace (it depends on the seed), keep the largest clean rectangle of at
   least 1 cm2 with ground already kept excluded. Reproduced the stored rectangle exactly on all four traces re-cut.
8. **Look**: render (28 planes along the normal), read, and draw (`code/make_qa.py`, `code/xsec_patch.py`).

Not in this recipe: the 28-plane renderer and the readers (project code and public weights), Arc/GPU reading, and the project's download meter. `fetch_mirror.py` and `trace_round.py` are cleaned copies that were
not re-run in this form.
