"""Resolve per-stain segmentation and object methods (--seg_method, --object_method).

Both options take a list of specs, each either ``STAIN=VALUE[,VALUE...]`` for one
stain or a bare ``VALUE[,VALUE...]`` for every stain not named. For each stain the
first value is its *primary* method; any others are comparisons.

Precedence, highest first: ``STAIN=`` on the command line, a bare value on the
command line, the built-in per-stain default, the built-in global default.

Imported by the Snakefile; kept free of Snakemake so it can be unit tested.
"""

import re
from dataclasses import dataclass, field

OBJECT_METHODS = ("watershed", "components", "none")

# Mask methods. ``p`` is the decimal separator, so names stay BIDS-safe.
_NUM = r"[0-9]+(?:p[0-9]+)?"
SEG_METHOD_PATTERNS = {
    "th": re.compile(r"^th[0-9]+$"),
    "otsu": re.compile(r"^otsu\+k[0-9]+i[0-9]+$"),
    "gmm": re.compile(rf"^gmm\+n[0-9]+k{_NUM}$"),
    "dog": re.compile(rf"^dog\+s{_NUM}r{_NUM}k{_NUM}$"),
    "lantern": re.compile(r"^lantern$"),
}


class SegPlanError(ValueError):
    """A --seg_method or --object_method spec that cannot be used."""


@dataclass
class StainPlan:
    """Methods for one stain; the first of each list is the primary."""

    stain: str
    seg_methods: list
    object_methods: list
    seg_source: str = "default"
    object_source: str = "default"

    @property
    def seg_method(self):
        return self.seg_methods[0]

    @property
    def object_method(self):
        return self.object_methods[0]

    @property
    def has_objects(self):
        return self.object_method != "none"

    def pairs(self, objects=True):
        """Every (seg_method, object_method) of this stain, primary first.

        With ``objects=False``, every seg method once (for mask-level outputs),
        paired with None.
        """
        if not objects:
            return [(m, None) for m in self.seg_methods]
        return [
            (m, o) for m in self.seg_methods for o in self.object_methods if o != "none"
        ]


@dataclass
class SegPlan:
    stains: dict = field(default_factory=dict)  # stain -> StainPlan

    def __iter__(self):
        return iter(self.stains.values())

    def __contains__(self, stain):
        return stain in self.stains

    def __getitem__(self, stain):
        return self.stains[stain]

    @property
    def stains_with_objects(self):
        return [s for s, p in self.stains.items() if p.has_objects]

    def stains_using(self, seg_method):
        """Stains that use ``seg_method`` (as primary or comparison)."""
        return [s for s, p in self.stains.items() if seg_method in p.seg_methods]

    def all_seg_methods(self):
        out = []
        for p in self:
            out += [m for m in p.seg_methods if m not in out]
        return out

    def describe(self):
        """One line per stain, for printing at start-up."""
        lines = []
        for p in self:
            seg = ",".join(p.seg_methods)
            obj = ",".join(p.object_methods)
            lines.append(
                f"  {p.stain}: seg_method {seg} ({p.seg_source}), "
                f"object_method {obj} ({p.object_source})"
            )
        return "\n".join(lines)

    def provenance(self):
        """List of dicts, one per (stain, seg method, object method)."""
        out = []
        for p in self:
            for m in p.seg_methods:
                for o in p.object_methods:
                    out.append(
                        {
                            "Stain": p.stain,
                            "SegMethod": m,
                            "ObjectMethod": o,
                            "Primary": m == p.seg_method and o == p.object_method,
                            "SegSource": p.seg_source,
                            "ObjectSource": p.object_source,
                        }
                    )
        return out


def parse_specs(specs, option):
    """Split specs into (bare values or None, {stain: values}).

    Each value list keeps its order and drops repeats. A stain named twice, or
    two bare specs, is an error, so precedence is never ambiguous.
    """
    bare = None
    by_stain = {}
    for spec in specs or []:
        spec = str(spec).strip()
        if not spec:
            continue
        if "=" in spec:
            stain, _, values = spec.partition("=")
            stain = stain.strip()
            if not stain:
                raise SegPlanError(f"{option} {spec}: missing stain before '='")
            if stain in by_stain:
                raise SegPlanError(f"{option}: stain {stain} is given more than once")
            by_stain[stain] = _split_values(values, option, spec)
        else:
            if bare is not None:
                raise SegPlanError(
                    f"{option}: more than one spec without a stain ({spec}); give "
                    "several methods as one comma-separated list instead"
                )
            bare = _split_values(spec, option, spec)
    return bare, by_stain


def _split_values(values, option, spec):
    out = []
    for v in values.split(","):
        v = v.strip()
        if not v:
            raise SegPlanError(f"{option} {spec}: empty method")
        if v not in out:
            out.append(v)
    return out


def check_seg_method(method):
    if method.endswith("+inst"):
        raise SegPlanError(
            f"--seg_method {method}: the '+inst' suffix is no longer supported. "
            f"Use --seg_method {method[: -len('+inst')]} with --object_method "
            "watershed (or STAIN=watershed)."
        )
    if not any(p.match(method) for p in SEG_METHOD_PATTERNS.values()):
        raise SegPlanError(
            f"--seg_method {method}: unknown method; expected one of th<value>, "
            "otsu+k<k>i<i>, gmm+n<n>k<k>, dog+s<s>r<r>k<k> or lantern"
        )


def check_object_method(method):
    if method not in OBJECT_METHODS:
        hint = " (connected components are now 'components')" if method == "cc" else ""
        hint = " (instances are now 'watershed')" if method == "inst" else hint
        raise SegPlanError(
            f"--object_method {method}: expected one of {', '.join(OBJECT_METHODS)}"
            f"{hint}"
        )


def _resolve_one(stain, cli_bare, cli_by_stain, default_by_stain, default):
    if stain in cli_by_stain:
        return list(cli_by_stain[stain]), "cli"
    if cli_bare is not None:
        return list(cli_bare), "cli"
    if stain in default_by_stain:
        value = default_by_stain[stain]
        return (list(value) if isinstance(value, list) else [value]), "default"
    return (list(default) if isinstance(default, list) else [default]), "default"


def torch_has_cuda():
    """Whether the installed torch is built with CUDA (the gpu pixi env).

    Reads torch/version.py rather than importing torch, which takes over a
    second and would slow every parse of the Snakefile. This reflects the
    environment, not the machine, so it is the same on a GPU-less submit node
    as on the GPU nodes its jobs run on.
    """
    import importlib.util
    from pathlib import Path

    spec = importlib.util.find_spec("torch")
    if spec is None or not spec.origin:
        return False
    try:
        text = Path(spec.origin).with_name("version.py").read_text()
    except OSError:
        return False
    return version_py_has_cuda(text)


def version_py_has_cuda(text):
    """Whether torch's version.py text names a CUDA version (``cuda = '12.6'``)."""
    pattern = r"^cuda\s*(?::[^=]*)?=\s*['\"][^'\"]+['\"]"
    return re.search(pattern, text, re.M) is not None


def apply_gpu_fallback(plan, gpu_available, gpu_fallbacks):
    """Replace GPU methods by their CPU fallbacks when there is no GPU.

    ``gpu_fallbacks`` maps each GPU seg method to its fallback. Methods that
    came from the built-in defaults are replaced (and reported); a GPU method
    asked for on the command line is an error, rather than a silent switch to
    another method or hours of CPU inference.

    Returns a list of notes describing the replacements.
    """
    notes = []
    if gpu_available:
        return notes
    for p in plan:
        gpu_here = [m for m in p.seg_methods if m in gpu_fallbacks]
        if not gpu_here:
            continue
        if p.seg_source == "cli":
            raise SegPlanError(
                f"--seg_method {p.stain}={','.join(p.seg_methods)}: "
                f"{', '.join(gpu_here)} needs a GPU, and this environment has no "
                "CUDA-enabled torch. Run in the gpu pixi environment (pixi run -e "
                "gpu ...), pass --gpu yes if the jobs run where a GPU is available, "
                f"or choose another method (e.g. {gpu_fallbacks[gpu_here[0]]})."
            )
        new = []
        for m in p.seg_methods:
            m = gpu_fallbacks.get(m, m)
            if m not in new:
                new.append(m)
        notes.append(
            f"{p.stain}: no GPU, so using {','.join(new)} instead of the default "
            f"{','.join(p.seg_methods)} (--gpu, gpu_seg_method_fallbacks)"
        )
        p.seg_methods = new
        p.seg_source = "default, no-GPU fallback"
    return notes


def resolve_seg_plan(
    stains,
    stains_for_seg,
    seg_specs,
    object_specs,
    seg_defaults,
    object_defaults,
):
    """Build the SegPlan for the stains of a dataset.

    Parameters
    ----------
    stains : list of str
        Stains in the data.
    stains_for_seg : list of str
        Stains to segment if present (--stains_for_seg). Stains named with
        ``STAIN=`` in either option are segmented too.
    seg_specs, object_specs : list of str or None
        --seg_method and --object_method as given (None if not given).
    seg_defaults, object_defaults : dict
        ``{"default": value, "stains": {stain: value}}`` built-in defaults.

    Returns
    -------
    (SegPlan, list of str)
        The plan, and notes worth printing (e.g. stains named but not in the data).
    """
    seg_bare, seg_by_stain = parse_specs(seg_specs, "--seg_method")
    obj_bare, obj_by_stain = parse_specs(object_specs, "--object_method")

    for m in (seg_bare or []) + [m for v in seg_by_stain.values() for m in v]:
        check_seg_method(m)
    for m in (obj_bare or []) + [m for v in obj_by_stain.values() for m in v]:
        check_object_method(m)

    notes = []
    named = list(dict.fromkeys(list(seg_by_stain) + list(obj_by_stain)))
    for s in named:
        if s not in stains:
            notes.append(
                f"stain {s} given in --seg_method/--object_method is not in the data {stains}"
            )

    to_segment = [s for s in stains if s in stains_for_seg or s in named]

    plan = SegPlan()
    for stain in to_segment:
        seg_methods, seg_source = _resolve_one(
            stain,
            seg_bare,
            seg_by_stain,
            seg_defaults.get("stains") or {},
            seg_defaults["default"],
        )
        obj_methods, obj_source = _resolve_one(
            stain,
            obj_bare,
            obj_by_stain,
            object_defaults.get("stains") or {},
            object_defaults["default"],
        )
        for m in seg_methods:
            check_seg_method(m)
        for m in obj_methods:
            check_object_method(m)
        if "none" in obj_methods and len(obj_methods) > 1:
            raise SegPlanError(
                f"--object_method for {stain}: 'none' cannot be combined with other "
                "object methods"
            )
        plan.stains[stain] = StainPlan(
            stain, seg_methods, obj_methods, seg_source, obj_source
        )
    return plan, notes
