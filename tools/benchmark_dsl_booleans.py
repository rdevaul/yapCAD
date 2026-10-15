"""Compare retained and eager DSL booleans without a private CAD fixture.

python tools/benchmark_dsl_booleans.py [--eager] [--command LOOP|BATCH]
The timer includes final mesh materialization, but excludes file IO/rendering.
"""

import argparse
import json
from pathlib import Path
import time

from yapcad.dsl import compile_and_run
from yapcad.dsl.runtime import builtins
from yapcad.boolean.manifold_engine import _to_manifold

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--eager", action="store_true")
parser.add_argument("--command", choices=("LOOP", "BATCH"), default="LOOP")
args = parser.parse_args()
if args.eager:
    # Reproduce the former eager path; this is benchmark-only instrumentation.
    builtins.try_boolean = lambda *args: None
source = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "performance"
    / "perforated_plate.dsl"
).read_text()
start = time.perf_counter()
result = compile_and_run(source, args.command, {})
elapsed = time.perf_counter() - start
if not result.success:
    raise RuntimeError(result.error_message)
body = _to_manifold(result.geometry, "benchmark")
print(
    json.dumps(
        dict(
            command=args.command,
            eager=args.eager,
            seconds=elapsed,
            volume=body.volume(),
            bbox=body.bounding_box(),
            triangles=body.num_tri(),
            status=body.status().name,
        )
    )
)
