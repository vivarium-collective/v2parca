#!/usr/bin/env python
"""
ParCa Pipeline as Process-Bigraph Steps
========================================

Runs the full vParCa pipeline — 9 Steps wired to a nested bigraph store
mirroring ``SimulationDataEcoli``'s structure — and pickles the final
store state.

Usage::

    # Fast run (reduced TF conditions via debug=True, ~30 min)
    python scripts/parca_bigraph.py --mode fast --cpus 4

    # Full production run (several hours)
    python scripts/parca_bigraph.py --mode full --cpus 4

    # Custom output directory and Km-fitting cache
    python scripts/parca_bigraph.py --mode fast --cpus 4 \\
        -o out/sim_data --cache-dir out/cache
"""

import argparse
import json
import os
import pickle
import re
import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from vparca.composite import build_parca_composite
from vparca.reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli


def main():
    parser = argparse.ArgumentParser(
        description="ParCa Pipeline as Process-Bigraph Steps",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--mode", choices=["fast", "full"], default="fast",
        help="fast: debug=True (reduced TF conditions, ~30 min); "
             "full: all conditions (several hours).",
    )
    parser.add_argument("-c", "--cpus", type=int, default=1,
                        help="Parallelism for Steps 4 and 5.")
    parser.add_argument("-o", "--outdir", type=str, default="out/sim_data",
                        help="Output directory for pickled state.")
    parser.add_argument("--cache-dir", type=str, default=None,
                        help="Cache directory for Km optimization "
                             "(default: <outdir>/cache).")
    parser.add_argument("--no-operons", action="store_true",
                        help="Disable operons in the KB.")
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    cache_dir = args.cache_dir or os.path.join(outdir, "cache")
    os.makedirs(cache_dir, exist_ok=True)

    print(f"\n{'=' * 60}\n  vParCa — {args.mode} mode\n{'=' * 60}")

    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] Loading raw_data (operons={not args.no_operons})")
    raw = KnowledgeBaseEcoli(
        operons_on=not args.no_operons,
        remove_rrna_operons=False, remove_rrff=False, stable_rrna=False,
    )
    print(f"    raw_data loaded in {time.time() - t0:.1f}s")

    t1 = time.time()
    print(f"\n[{time.strftime('%H:%M:%S')}] Running ParCa pipeline ...")
    # Capture per-step timings by tee-ing the Step completion lines the
    # Steps already emit ("  Step N (name) completed in X.Xs").
    import io
    import contextlib
    captured = io.StringIO()
    class _Tee:
        def __init__(self, *streams): self.streams = streams
        def write(self, s):
            for st in self.streams: st.write(s)
        def flush(self):
            for st in self.streams: st.flush()
    with contextlib.redirect_stdout(_Tee(sys.stdout, captured)):
        composite = build_parca_composite(
            raw,
            debug=(args.mode == 'fast'),
            cpus=args.cpus,
            cache_dir=cache_dir,
        )
    print(f"\n[{time.strftime('%H:%M:%S')}] Pipeline completed in "
          f"{time.time() - t1:.1f}s")

    # Parse per-step timings from the captured stdout.
    step_times = {}
    for m in re.finditer(r'Step (\d) .*? completed in ([0-9.]+)s',
                         captured.getvalue()):
        step_times[f'step_{m.group(1)}'] = float(m.group(2))

    runtimes_path = os.path.join(outdir, 'runtimes.json')
    with open(runtimes_path, 'w') as f:
        json.dump(step_times, f, indent=2, sort_keys=True)
    print(f"    runtimes: {runtimes_path}")

    # Pickle the full store state — subsystem objects + top-level dicts +
    # cell_specs entries.  Strip the tick leaves (sequencing tokens, not
    # useful downstream) and the internal bigraph-schema bookkeeping keys.
    state = {k: v for k, v in composite.state.items()
             if not k.startswith('tick_')
             and not k.startswith('_')
             and k not in {'global_time', 'initialize', 'input_adjustments',
                           'basal_specs', 'tf_condition_specs', 'fit_condition',
                           'promoter_binding', 'adjust_promoters',
                           'set_conditions', 'final_adjustments'}}

    out_path = os.path.join(outdir, "parca_state.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
    size_mb = os.path.getsize(out_path) / (1024 * 1024)

    total = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"Total time:  {total:.1f}s ({total / 60:.1f} min)")
    print(f"Output:      {out_path} ({size_mb:.1f} MB)")
    print(f"Store keys:  {sorted(state.keys())}")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()
