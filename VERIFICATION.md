# Release verification

The complete primary workflow was rerun successfully on 2026-09-29:

```bash
python3 verify.py --jobs 8
```

Final result: `FULL CERTIFICATE REPLAY PASSED: all.`

## Environment

- macOS, arm64
- Python 3.12.14; NumPy 2.3.5
- GCC 15.2.0; MPFR 4.2.2; GMP 6.3.0
- C++17, `-O3 -fno-fast-math -ffp-contract=off -pthread`

## Results

| Component | Result |
| --- | --- |
| Unrestricted exact subdivision | 167,505 leaves; complete binary partition |
| Unrestricted binary export | Byte-for-byte agreement |
| Unrestricted nonlocal interval checks | 166,505 boxes; zero failures |
| Unrestricted local and endpoint checks | 512 Hessian boxes; all boundary, value, and localization checks passed |
| MHR exact interior subdivision | 64,175 leaves; prefix-free; exact Kraft sum 1 |
| MHR interior interval replay | 64,175 leaves; zero failures |
| MHR local checks | Strict convexity and both refined residual signs passed |
| MHR exterior replay | 47,273,387 nodes; 23,636,694 leaves; zero failures |

The release also checked that removing a leaf is rejected by the MHR
geometry audit and that the runner refuses disabled Python assertions or
a nonpositive worker count. All five C++ programs compiled without warnings
in the environment above.

The finite certificate data and mathematical verification logic are
unchanged from the manuscript's accompanying archive. The repository adds
a portable runner and removes the archive-integrity wrapper and the
geometry audit's informational digest output. Three source comments were
updated to refer to the companion paper.

The optional mpmath workflow was not rerun for this packaging revision.
Successful replay checks the stated finite certificate against the
acceptance inequalities; the analytic justification is in the paper.
