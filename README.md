# Sharp Welfare Certificates

Reproducible computer-assisted proof data for the paper

> **Sharp Second-Best Welfare in Bilateral and Matching Markets**

This repository contains the finite proof inputs and the verification programs
that certify the two sharp welfare constants of the paper. It is the
standalone, self-contained release of the certificate archive described in
Appendix D of the manuscript; no other files from the submission are needed
to run it.

## Results being certified

The paper determines how much total welfare a market must lose because values
are private. Under independent nonnegative types, Bayesian incentive
compatibility (BIC), interim individual rationality (IR), and ex ante weak
budget balance (WBB), the sharp ratio of second-best to first-best welfare is:

| Buyer priors | Certified ratio $\rho$ |
| --- | --- |
| Arbitrary independent priors | `.8882516903` < $\rho_{\mathrm U}$ < `.8882516904` |
| MHR buyers; arbitrary sellers | `.911389368124` < $\rho_{\mathrm M}$ < `.911389368127` |

Exact interval brackets, as printed in Theorem 1.1 of the paper:

| Constant | Certified enclosure |
| --- | --- |
| $\rho_{\mathrm U}$ | `0.88825169032566246988432630129334211037 < rho_U < 0.88825169032566246988432630129334211038` |
| $\rho_{\mathrm M}$ | `0.911389368124 < rho_M < 0.911389368127` |

Both constants are already attained in bilateral trade. The unrestricted-prior
guarantee holds for every downward-closed family of feasible matchings; the MHR
guarantee holds for ordinary matching markets, where every matching of the
compatibility graph is feasible. Buyer values are scalar, and matches are
evaluated with unrestricted signed transfers. The same guarantees and sharp
constants hold under pointwise strong budget balance.

## What this repository does and does not do

These programs **certify the two bilateral welfare constants** by replaying
finite interval certificates against the acceptance inequalities derived in
Appendices B and C of the paper.

They do **not** enumerate market instances, and they do not re-derive the
matching-market results. The paper's analytic lifting arguments — a welfare
endowment charge for unrestricted priors, and an MHR-preserving rectangle
packing argument for MHR buyers — transfer the bilateral inequalities to the
stated matching markets without loss. Those arguments are proved in the
manuscript; the finite computation here supplies their numerical foundation.

The proof is computer assisted, but in a restricted and checkable sense. It
does not rely on an optimizer's termination, on floating-point sign guesses,
or on a formal proof assistant. All analytic reductions, acceptance
inequalities, and remainder bounds appear in the manuscript; this repository
supplies their finite inputs and implementations.

## Run the complete verification

**Requirements.** Python 3.10 or later, NumPy, a C++17 compiler with the
MPFR/GMP development libraries, and a little-endian machine. Python assertions
must remain enabled (do not pass `-O` and do not set `PYTHONOPTIMIZE`); the
runner enforces this.

On Ubuntu/Debian:

```bash
sudo apt-get install g++ libmpfr-dev libgmp-dev python3-venv
```

On macOS with Homebrew:

```bash
brew install gcc mpfr gmp
```

Then, from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 verify.py --jobs 8
```

The runner compiles the C++ verifiers, expands the compressed exterior proof
tree into `build/`, and runs every check. It stops at the first failure and
writes each program's output to `logs/`. A successful complete run ends with:

```text
FULL CERTIFICATE REPLAY PASSED: all.
```

The full replay can take several minutes depending on the machine. The worker
count affects only the MHR checks. To run a single certificate:

```bash
python3 verify.py --suite unrestricted
python3 verify.py --suite mhr --jobs 8
```

`CXX` selects the compiler; `CPPFLAGS`, `CXXFLAGS`, and `LDFLAGS` add build
options. On macOS the runner detects Homebrew library paths and prefers an
installed Homebrew GCC. It always disables fast-math and floating-point
contraction for the interval calculations, and these flags are appended after
any user flags so they cannot be overridden.

For a quick structural check that does no compilation and no interval
evaluation:

```bash
python3 verify.py --check-only
```

This checks subdivision structure and binary export consistency only. Its
success message says so explicitly — the interval proofs have not run.

## What is checked

**Unrestricted priors.** `verify_geometry.py` reconstructs the exact binary
partition of 167,505 leaves. `verify_export.py` checks that `cover_boxes.bin`
is a byte-for-byte export of `cover.npz`. The MPFR verifier then checks all
166,505 nonlocal boxes. A second MPFR program checks the exterior
inequalities, the 512 local Hessian boxes, and the final value, boundary, and
localization bounds.

**MHR buyers.** `audit_interior.py` checks the complete, prefix-free interior
subdivision of 64,175 leaves using exact integer arithmetic. The replay
program verifies every interior leaf, and the exterior program traverses all
47,273,387 tree nodes and checks the exterior inequalities. The local program
verifies strict convexity and the refined residual signs.

All of these components must pass. The local program's own success message,
`PASS LOCAL CHECKS ONLY`, refers to the local component alone and is **not** a
global certificate — the global claim also requires both covers.

The input files propose finite subdivisions; their labels are not accepted as
mathematical signs. The verifiers recompute enclosing inequalities rather than
trusting the labels. The coverage checks, the binary export comparison, and
the interval acceptance tests are all retained in this version.

## Repository layout

```text
verify.py                  Complete verification entry point
requirements.txt           Dependencies for the primary MPFR workflow
requirements-optional.txt  Optional Python checks and cover generation
unrestricted/              Power-law certificate, data, and verifiers
mhr/                       MHR certificate, data, and verifiers
VERIFICATION.md            Results of the release verification run
```

The required finite inputs are committed directly, so no preliminary archive
extraction is needed. `build/`, `logs/`, and Python caches are generated
locally and ignored by Git. The compressed MHR exterior tree is expanded
automatically, leaving its tracked copy unchanged.

## Optional Python implementation

The unrestricted certificate also ships an mpmath interval implementation and
auxiliary consistency checks. These are useful for inspection and for
cross-checking the C++ programs, but they are slower and are **not** a
replacement for the interval proof. From the repository root:

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

The `0` in the cover command requests every nonlocal box.

`generate_cover.py` can regenerate the unrestricted subdivision using SciPy.
It overwrites `cover.npz`, and it is unnecessary for verifying the supplied
data. If you do regenerate, re-export the binary with
`python3 export_cover_binary.py` and rerun the complete verification.

## Verification status

See [`VERIFICATION.md`](VERIFICATION.md) for the release test run: the
complete primary workflow was rerun successfully on 2026-09-29, ending in
`FULL CERTIFICATE REPLAY PASSED: all.` That run also confirmed that removing a
leaf is rejected by the MHR geometry audit, and that the runner refuses
disabled Python assertions or a nonpositive worker count.

## Citation

If you use these certificates, please cite the paper. The manuscript is
currently under submission; citation details will be added here once they are
available.

## License

No license has been chosen yet, so all rights are reserved by default. If you
would like to reuse this material, please contact the author, or watch for a
license to be added.
