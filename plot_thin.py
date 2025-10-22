#!/usr/bin/env python3

import cloudpickle
import sys
from dask.diagnostics import visualize, Profiler, ResourceProfiler, CacheProfiler
import bokeh.plotting as bp

# -----------------------------
# Load profiling data
# -----------------------------
if len(sys.argv) < 2:
    print("Usage: python profile_to_mov.py <profile.pkl>")
    sys.exit(1)

pkl = sys.argv[1]

with open(pkl, "rb") as f:
    prof, rprof, cprof = cloudpickle.load(f)


# 🎨 Visualize thinned data
pt = visualize(prof, filename=f"/nfs/khan/trainees/akhan488/profile_task.html",save=False,show=False)
bp.output_file(f"/nfs/khan/trainees/akhan488/profile_tasks.html")
bp.save(pt)
pr = visualize(rprof,save=False,show=False)
bp.output_file(f"/nfs/khan/trainees/akhan488/profile_resources.html")
bp.save(pr)

ptr = visualize([prof,rprof],save=False,show=False)
bp.output_file(f"/nfs/khan/trainees/akhan488/profile_tasks+resources.html")
bp.save(ptr)



print(f"✅ Visualization saved")

