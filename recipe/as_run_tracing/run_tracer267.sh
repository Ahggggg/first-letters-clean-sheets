#!/usr/bin/env bash
# T-267: one detached tracer container (no window).  usage: run_tracer267.sh NAME VOLUME GENERATIONS X Y Z [CPUS]
# VOLUME = docker named volume holding a copy of the local mirror (pred.zarr/, ngrid/) on the VM's ext4 disk (a Windows bind mount made the tracer I/O-bound:
# 25 generations in 28 min at 33 pct of one core, process in state D).  Mirror volume -> /data (read-only); output -> data/t267-trace/runs/NAME/out (segment dir "seg").
# Params = the official PHerc0826 auto_grown_20260902140841364 vc_gsfs_params with generations set here and cache_size 3e9, cache_root inside the container.
# Container name t267-NAME; poll with docker logs / docker ps; remove with docker rm after the result is read.  Entry point overridden (no AWS credentials needed).
set -euo pipefail
export PATH="<local path>:$PATH"
export MSYS_NO_PATHCONV=1
NAME=$1; VOL=$2; GEN=$3; X=$4; Y=$5; Z=$6; CPUS=${7:-1.5}
ROOT=<local path>
RUN=$ROOT/runs/$NAME
mkdir -p "$RUN/out"
cat > "$RUN/params.json" <<EOF
{"cache_root": "/tmp/cache", "cache_size": 3000000000, "generations": $GEN, "min_area_cm": 0.3, "mode": "seed", "normal_grid_path": "/data/ngrid", "thread_limit": 1, "use_cuda": false, "voxelsize": 9.362}
EOF
W() { cygpath -w "$1"; }
docker run -d --name "t267-$NAME" --cpus "$CPUS" --memory 8g \
  -v "$VOL:/data:ro" -v "$(W $RUN/params.json):/params.json:ro" -v "$(W $RUN/out):/out" \
  --entrypoint vc_grow_seg_from_seed ghcr.io/scrollprize/villa/volume-cartographer:main \
  -v /data/pred.zarr -t /out/seg -p /params.json --seed "$X" "$Y" "$Z" --segment-name seg
