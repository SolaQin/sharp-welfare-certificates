#!/usr/bin/env python3
"""Build and replay the finite certificates for both sharp welfare constants."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build"
LOGS = ROOT / "logs"


def run(name, command, cwd=ROOT):
    print(f"Running {name} ...", flush=True)
    log = LOGS / f"{name}.log"
    with log.open("w") as output:
        result = subprocess.run(command, cwd=cwd, stdout=output,
                                stderr=subprocess.STDOUT)
    if result.returncode:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-20:])
        raise RuntimeError(f"{name} failed (exit {result.returncode}).\n{tail}\n"
                           f"Full output: {log}")
    print(f"Completed {name}; output in logs/{log.name}", flush=True)


def compiler_options():
    default = "c++"
    if sys.platform == "darwin":
        candidates = []
        for prefix in (Path("/opt/homebrew/bin"), Path("/usr/local/bin")):
            candidates.extend(p for p in prefix.glob("g++-*")
                              if p.is_file() and p.name[4:].isdigit())
        if candidates:
            default = str(max(candidates, key=lambda p: int(p.name[4:])))
    compiler = shlex.split(os.environ.get("CXX", default))
    if not compiler or not shutil.which(compiler[0]):
        raise RuntimeError("C++ compiler not found. Install a C++17 compiler or set CXX.")
    flags = shlex.split(os.environ.get("CPPFLAGS", ""))
    flags += shlex.split(os.environ.get("CXXFLAGS", "-O3"))
    libraries = shlex.split(os.environ.get("LDFLAGS", ""))
    if sys.platform == "darwin":
        for prefix in (Path("/opt/homebrew"), Path("/usr/local")):
            for package in ("mpfr", "gmp"):
                location = prefix / "opt" / package
                if (location / "include").is_dir():
                    flags.append(f"-I{location / 'include'}")
                    libraries.append(f"-L{location / 'lib'}")
    # These flags come last so user optimization flags cannot enable contraction
    # or fast arithmetic in the interval calculations.
    flags += ["-std=c++17", "-fno-fast-math", "-ffp-contract=off", "-pthread"]
    libraries += ["-lmpfr", "-lgmp"]
    return compiler, flags, libraries


def compile_one(directory, source, name, options):
    compiler, flags, libraries = options
    target = BUILD / name
    run(name + "_compile", [*compiler, *flags, source, *libraries,
                            "-o", str(target)], directory)
    return str(target)


def unrestricted(check_only, options):
    directory = ROOT / "unrestricted"
    run("unrestricted_geometry", [sys.executable, "verify_geometry.py"], directory)
    run("unrestricted_export", [sys.executable, "verify_export.py"], directory)
    if check_only:
        return
    cover = compile_one(directory, "verify_cover_mpfr.cpp", "unrestricted_cover", options)
    local = compile_one(directory, "verify_local_point_mpfr.cpp", "unrestricted_local", options)
    run("unrestricted_cover", [cover, "cover_boxes.bin", "0"], directory)
    run("unrestricted_local", [local], directory)


def mhr(check_only, options, jobs):
    directory = ROOT / "mhr"
    run("mhr_geometry", [sys.executable, "audit_interior.py", "interior.leaves"], directory)
    if check_only:
        return
    replay = compile_one(directory, "replay_residual_cover.cpp", "mhr_replay", options)
    local = compile_one(directory, "verify_candidate_neighborhood.cpp", "mhr_local", options)
    exterior = compile_one(directory, "verify_exterior_cover.cpp", "mhr_exterior", options)
    run("mhr_local", [local], directory)
    run("mhr_interior", [replay, "interior.leaves", str(jobs)], directory)
    data = BUILD / "exterior_cover.bin"
    with gzip.open(directory / "exterior_cover.bin.gz", "rb") as source:
        with data.open("wb") as target:
            shutil.copyfileobj(source, target)
    run("mhr_exterior", [exterior, str(data), ".05", str(jobs)], directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("all", "unrestricted", "mhr"), default="all")
    parser.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                        help="worker count for the MHR replay (default: up to 8)")
    parser.add_argument("--check-only", action="store_true",
                        help="check subdivision/export structure only; do not certify the ratios")
    args = parser.parse_args()
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        parser.error("Python assertions must remain enabled; unset PYTHONOPTIMIZE and omit -O.")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if sys.byteorder != "little":
        parser.error("The supplied binary readers require a little-endian machine.")
    if args.suite in ("all", "unrestricted"):
        try:
            import numpy  # noqa: F401
        except ImportError:
            parser.error("NumPy is required: python3 -m pip install -r requirements.txt")
    options = None if args.check_only else compiler_options()
    BUILD.mkdir(exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    if args.suite in ("all", "unrestricted"):
        unrestricted(args.check_only, options)
    if args.suite in ("all", "mhr"):
        mhr(args.check_only, options, args.jobs)
    if args.check_only:
        print("STRUCTURAL CHECKS PASSED. Interval proofs have not been run.")
    else:
        print(f"FULL CERTIFICATE REPLAY PASSED: {args.suite}.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
