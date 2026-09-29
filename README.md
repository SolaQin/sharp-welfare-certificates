# Sharp Welfare Certificates

Verification code and finite proof data for **Sharp Second-Best Welfare in
Bilateral and Matching Markets**. The analytic reductions and acceptance
inequalities are given in Appendices B and C of the paper.

The programs certify the two bilateral welfare constants. The paper's
analytic lifting arguments transfer these constants to the stated matching
markets; the programs do not enumerate market instances.

| Buyer priors | Certified welfare ratio |
| --- | --- |
| Arbitrary independent priors | `0.88825169032566246988432630129334211037 < rho_U < 0.88825169032566246988432630129334211038` |
| MHR buyers; arbitrary independent sellers | `0.911389368124 < rho_M < 0.911389368127` |

## Run the complete verification

Requirements: Python 3.10 or later, NumPy, a C++17 compiler, and the MPFR/GMP
development libraries. The binary readers use little-endian data. Python
assertions must remain enabled.

On Ubuntu/Debian, install the native dependencies with:

```bash
sudo apt-get install g++ libmpfr-dev libgmp-dev python3-venv
```

On macOS with Homebrew:

```bash
brew install gcc mpfr gmp
```

From this repository's root directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 verify.py --jobs 8
```

The runner compiles the verifiers, expands the exterior proof tree into
`build/`, and runs all checks. It stops on the first failure and writes each
program's output to `logs/`. A successful complete run ends with:

```text
FULL CERTIFICATE REPLAY PASSED: all.
```

The full replay can take several minutes, depending on the machine. The
worker count affects the MHR checks only. To run one certificate:

```bash
python3 verify.py --suite unrestricted
python3 verify.py --suite mhr --jobs 8
```

`CXX` selects the compiler; `CPPFLAGS`, `CXXFLAGS`, and `LDFLAGS` add build
options. On macOS the runner detects Homebrew library paths and prefers an
installed Homebrew GCC. It always disables fast-math and floating-point
contraction for the interval calculations.

## What is checked

**Unrestricted priors.** `verify_geometry.py` reconstructs the exact binary
partition of 167,505 leaves. `verify_export.py` checks that `cover_boxes.bin`
is a byte-for-byte export of `cover.npz`. The MPFR verifier checks all 166,505
nonlocal boxes. The second MPFR program checks the exterior inequalities,
512 local Hessian boxes, and the final value, boundary, and localization
bounds.

**MHR buyers.** `audit_interior.py` checks the complete, prefix-free interior
subdivision of 64,175 leaves using exact integer arithmetic. The replay
program verifies every interior leaf, and the exterior program traverses
all 47,273,387 tree nodes and checks the exterior inequalities. The local
program verifies strict convexity and the refined residual signs. All these
components must pass: `PASS LOCAL CHECKS ONLY` is not a global certificate.

The input files propose finite subdivisions; their labels are not accepted
as mathematical signs. The verifiers recompute enclosing inequalities. The
coverage checks, binary export comparison, and interval acceptance tests
are all retained in this version.

For a quick structural check without compilation or interval evaluation:

```bash
python3 verify.py --check-only
```

This checks subdivision structure and binary export consistency only. Its
success message explicitly states that the interval proofs have not run.

## Repository layout

```text
verify.py                  Complete verification entry point
requirements.txt           Dependencies for the primary MPFR workflow
requirements-optional.txt  Optional Python checks and cover generation
unrestricted/              Power-law certificate, data, and verifiers
mhr/                       MHR certificate, data, and verifiers
VERIFICATION.md            Results of the release verification run
```

The repository contains the required finite inputs directly; no preliminary
archive extraction is needed. `build/`, `logs/`, and Python caches are
generated locally and ignored by Git. The compressed MHR exterior tree is
expanded automatically, leaving its tracked copy unchanged.

## Optional Python implementation

The unrestricted certificate also includes an mpmath interval implementation
and auxiliary consistency checks. Starting at the repository root, install
their dependencies and then run the checks in `unrestricted/`:

```bash
python3 -m pip install -r requirements-optional.txt
cd unrestricted
python3 verify_geometry.py
python3 verify_boundaries.py
python3 verify_cover.py cover.npz 0 8
python3 verify_local.py
python3 verify_point.py
python3 verify_export.py
python3 check_consistency.py
```

The `0` in the cover command requests every nonlocal box. The auxiliary
consistency checks are not a replacement for the interval proof.
`generate_cover.py` can regenerate the unrestricted subdivision using
SciPy; it overwrites `cover.npz` and is unnecessary for verification of the
supplied data. After regeneration, export the new binary data with
`python3 export_cover_binary.py` and rerun the complete verification.

## Upload to GitHub

Use this directory's contents as the repository root. The included data
files are ready to upload with the source; no compiled binaries, local
environments, or generated logs are required.
