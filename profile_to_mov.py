#!/usr/bin/env python3
"""
Render Dask Profiler (.pkl via cloudpickle) into an animated GIF and MP4.
Each frame shows a sliding time window of task execution.

Usage:
    python profile_to_mov.py profile.pkl
"""

import sys, os, tempfile, subprocess
import numpy as np
import cloudpickle
from dask.diagnostics import visualize
from PIL import Image
import imageio.v2 as imageio

# -----------------------------
# Parameters
# -----------------------------
window_size = 5.0   # seconds per visualization window
step_size   = 2.0   # slide step in seconds
fps         = 2     # frames per second for GIF/video
outfile_gif = "profile_animation.gif"
outfile_mp4 = "profile_animation.mp4"
chromium_bin = "chromium-browser"   # or "google-chrome" on some systems

# -----------------------------
# Load profiling data
# -----------------------------
if len(sys.argv) < 2:
    print("Usage: python profile_to_mov.py <profile.pkl>")
    sys.exit(1)

pkl = sys.argv[1]
with open(pkl, "rb") as f:
    prof, rprof, *rest = cloudpickle.load(f)

print(f"Loaded profiler with {len(prof.results)} task records")
print(f"Loaded resource profiler with {len(rprof.results)} samples")

# -----------------------------
# Helpers
# -----------------------------
def slice_prof(prof, t_start, t_end):
    """Return a new profiler with results in [t_start, t_end)."""
    sliced = prof.__class__()
    sliced.results = [
        p for p in prof.results
        if p.end_time >= t_start and p.start_time < t_end
    ]
    return sliced

def slice_rprof(rprof, t_start, t_end):
    """Return a new ResourceProfiler containing only samples in [t_start, t_end)."""
    sliced = rprof.__class__(dt=rprof._dt)
    sliced.results = [r for r in rprof.results if t_start <= r.time < t_end]
    return sliced



def render_html_to_image(html_path):
    """Render a Bokeh HTML visualization to PNG using headless Chromium."""
    png_path = tempfile.mktemp(suffix=".png")
    cmd = [
        chromium_bin, "--headless", "--disable-gpu",
        f"--screenshot={png_path}", "--window-size=1400,900", html_path
    ]
    try:
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return Image.open(png_path)
    except Exception as e:
        print(f"⚠️ Could not render HTML to image: {e}")
        return None

# -----------------------------
# Determine overall time range
# -----------------------------
t0 = min(p.start_time for p in prof.results)
t1 = max(p.end_time for p in prof.results)
print(f"Profiling duration: {t1 - t0:.2f} seconds")

# -----------------------------
# Generate frames
# -----------------------------
frames = []

with tempfile.TemporaryDirectory() as tmp:
    for t in np.arange(t0, t1, step_size):
        t_start, t_end = t, t + window_size
        prof_sub  = slice_prof(prof,  t_start, t_end)
        rprof_sub = slice_rprof(rprof, t_start, t_end)
        html_file = os.path.join(tmp, f"frame_{t_start:.2f}.html")

        # Skip windows with no data
        if not prof_sub.results and not rprof_sub.results:
            print(f"⏭️  Skipping empty window {t_start:.1f}-{t_end:.1f}s (no tasks or samples)")
            continue

        visualize([prof_sub, rprof_sub], filename=html_file, show=False)


        img = render_html_to_image(html_file)
        if img:
            frames.append(img.convert("RGB"))
            print(f"🖼️ Frame {len(frames):03d}: {t_start:.1f}-{t_end:.1f}s "
                  f"({len(prof_sub.results)} tasks, {len(rprof_sub.results)} samples)")
        else:
            print(f"⚠️ Skipped frame {t_start:.1f}-{t_end:.1f}s")

# -----------------------------
# Save animation
# -----------------------------
if not frames:
    print("❌ No frames captured. Make sure Chromium is installed and in PATH.")
    sys.exit(1)

print(f"Saving {len(frames)} frames to GIF ({outfile_gif}) and MP4 ({outfile_mp4})...")
imageio.mimsave(outfile_gif, frames, fps=fps)
imageio.mimsave(outfile_mp4, frames, fps=fps)
print("✅ Done.")


