"""Parse the MATLAB robustness ``Results.txt`` table of GDP moments.

The MATLAB replication code writes one block per reallocation mode and step
size, with rows of the form::

    (sigma, epsilon, theta) * mean * std * skew * exkurt * amplification

This module extracts the rows for the paper benchmark elasticities
(sigma=0.9, epsilon=0.5, theta=0.001) and prints them. Optionally, a JSON file
with Julia moments can be supplied via ``--julia-results`` to print a
side-by-side comparison; no Julia numbers are hard-coded here.

Expected ``--julia-results`` JSON layout (all fields optional)::

    {
      "No Reallocation|1":   {"mean": -0.0035, "std": 0.0112, "skew": -0.15, "exkurt": 0.12},
      "Full Reallocation|4": {"mean": -0.0097, "std": 0.0248}
    }

Run ``python bf_replication/parse_results.py --help`` for usage.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_RESULTS_PATH = (
    REPO_ROOT
    / "Replication Files"
    / "GDP Simulatin -- 88 Sector"
    / "Robustness"
    / "Results.txt"
)

#: Paper benchmark elasticities (sigma, epsilon, theta).
BENCHMARK_SIGMA = 0.9
BENCHMARK_EPSILON = 0.5
BENCHMARK_THETA = 0.001

#: (mode, step) rows that must be present for the script to succeed.
REQUIRED_ROWS: Tuple[Tuple[str, int], ...] = (
    ("No Reallocation", 1),
    ("No Reallocation", 4),
    ("Full Reallocation", 1),
    ("Full Reallocation", 4),
)

#: Strict numeric token: optional sign, mandatory digits in the mantissa, and a
#: well-formed optional exponent. This rejects malformed combinations such as
#: ``e+-`` (which would otherwise slip through ``[\d.eE+-]+``).
_NUM_RE = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

_ROW_RE = re.compile(
    r"\(" + "(" + _NUM_RE + ")" + r",\s*" + "(" + _NUM_RE + ")" + r",\s*" + "(" + _NUM_RE + ")" + r"\)"
    + r"\s*\*\s*" + "(" + _NUM_RE + ")"
    + r"\s*\*\s*" + "(" + _NUM_RE + ")"
    + r"\s*\*\s*" + "(" + _NUM_RE + ")"
    + r"\s*\*\s*" + "(" + _NUM_RE + ")"
    + r"\s*\*\s*" + "(" + _NUM_RE + ")"
    + r"\s*$"
)


class Moments(NamedTuple):
    """Moments of the simulated log-GDP distribution for one parameter row."""

    mean: float
    std: float
    skew: float
    exkurt: float
    amplification: float


def read_lines(path: Path) -> List[str]:
    """Return the lines of ``path``; raises ``FileNotFoundError`` if absent."""
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def find_moments(
    lines: Sequence[str],
    mode: str,
    step: int,
    sigma: float,
    epsilon: float,
    theta: float,
) -> Optional[Moments]:
    """Find the row matching the given mode, step size and elasticities.

    ``lines`` is any sequence of raw text lines from ``Results.txt``. Returns
    ``None`` when no matching row exists.
    """
    curr_mode = ""
    curr_step = ""
    for raw_line in lines:
        line = raw_line.strip()
        if "Reallocation" in line:
            curr_mode = line
            continue
        if line.startswith("Step size"):
            curr_step = line.split()[-1]
            continue
        if curr_mode != mode or curr_step != str(step):
            continue
        match = _ROW_RE.match(line)
        if not match:
            continue
        try:
            values = [float(match.group(i)) for i in range(1, 9)]
        except ValueError:
            # A matched row whose tokens do not parse should never happen with
            # the strict numeric regex, but skip it defensively instead of
            # raising an uncaught traceback.
            continue
        s, e, t = values[0], values[1], values[2]
        if abs(s - sigma) < 0.01 and abs(e - epsilon) < 0.01 and abs(t - theta) < 0.001:
            return Moments(*values[3:8])
    return None


#: Moment fields recognized from a Julia results entry.
_JULIA_MOMENT_FIELDS = ("mean", "std", "skew", "exkurt")


def load_julia_results(path: Path) -> Dict[str, Dict[str, float]]:
    """Load and validate a JSON object of Julia moments.

    The top-level value must be an object keyed by ``"<mode>|<step>"``; each
    entry must itself be an object/dict, and any recognized moment value must
    be a real number (``int``/``float``, but not ``bool``). Malformed input
    raises ``ValueError`` so the caller can emit a concise CLI error.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("top-level JSON must be an object/dict")
    for key, entry in data.items():
        if not isinstance(entry, dict):
            raise ValueError(f"entry {key!r} must be an object/dict")
        for field in _JULIA_MOMENT_FIELDS:
            if field not in entry:
                continue
            value = entry[field]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(
                    f"entry {key!r} field {field!r} must be numeric, "
                    f"got {type(value).__name__}"
                )
    return data


def format_matlab(moments: Moments) -> str:
    """Format MATLAB moments as a single report line."""
    return (
        f"mean={moments.mean * 100:.3f}%  std={moments.std * 100:.3f}%  "
        f"skew={moments.skew:.4f}  exkurt={moments.exkurt:.4f}  "
        f"ampl={moments.amplification:.4f}"
    )


def report_row(
    lines: Sequence[str],
    mode: str,
    step: int,
    julia: Optional[Dict[str, Dict[str, float]]],
) -> Optional[Moments]:
    """Print the MATLAB row (and any Julia comparison) for ``mode``/``step``."""
    header = (
        f"Paper benchmark (sigma={BENCHMARK_SIGMA}, epsilon={BENCHMARK_EPSILON}, "
        f"theta={BENCHMARK_THETA}), {mode}, step={step}:"
    )
    print(header)
    moments = find_moments(
        lines, mode, step, BENCHMARK_SIGMA, BENCHMARK_EPSILON, BENCHMARK_THETA
    )
    if moments is None:
        print("  MATLAB: row not found")
        return None

    print(f"  MATLAB: {format_matlab(moments)}")
    if julia is None:
        return moments

    entry = julia.get(f"{mode}|{step}")
    if not entry:
        print("  Julia:  no entry in --julia-results for this row")
        return moments

    fields = (("mean", 100.0, "pp"), ("std", 100.0, "pp"), ("skew", 1.0, ""), ("exkurt", 1.0, ""))
    julia_parts: List[str] = []
    diff_parts: List[str] = []
    for name, scale, unit in fields:
        value = entry.get(name)
        if value is None:
            continue
        matlab_value = getattr(moments, name)
        suffix = "%" if unit == "pp" else ""
        julia_parts.append(f"{name}={value * scale:.3f}{suffix}")
        diff_parts.append(f"{name}={abs(matlab_value - value) * scale:.4f}{unit}")
    if julia_parts:
        print(f"  Julia:  {'  '.join(julia_parts)}")
        print(f"  DIFF:   {'  '.join(diff_parts)}")
    else:
        print("  Julia:  entry present but contains no known moment fields")
    return moments


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse the MATLAB robustness Results.txt table of GDP moments."
    )
    parser.add_argument(
        "results",
        nargs="?",
        type=Path,
        default=DEFAULT_RESULTS_PATH,
        help="path to Results.txt (default: the in-repository robustness table)",
    )
    parser.add_argument(
        "--julia-results",
        type=Path,
        default=None,
        metavar="JSON",
        help=(
            "optional JSON file with Julia moments keyed by '<mode>|<step>'; "
            "without it only MATLAB values are printed"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point: returns 0 on success, 1 when input or rows are missing."""
    args = parse_args(argv)

    try:
        lines = read_lines(args.results)
    except OSError as exc:
        print(f"error: cannot read {args.results}: {exc}", file=sys.stderr)
        return 1

    julia: Optional[Dict[str, Dict[str, float]]] = None
    if args.julia_results is not None:
        try:
            julia = load_julia_results(args.julia_results)
        except (OSError, ValueError) as exc:
            print(
                f"error: cannot read Julia results {args.julia_results}: {exc}",
                file=sys.stderr,
            )
            return 1

    print(f"Source: {args.results}")
    if julia is None:
        print(
            "Julia comparison skipped: pass --julia-results with an external "
            "JSON artifact to compare against MATLAB."
        )
    print()

    missing: List[str] = []
    for mode, step in REQUIRED_ROWS:
        if report_row(lines, mode, step, julia) is None:
            missing.append(f"{mode}, step={step}")
        print()

    if missing:
        print(
            "error: required rows missing from "
            f"{args.results}: {'; '.join(missing)}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
