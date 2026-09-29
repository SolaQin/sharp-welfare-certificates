"""Read-only, exact audit of the finite binary subdivision cover."""
from pathlib import Path
from collections import Counter
import sys

path = Path(sys.argv[1] if len(sys.argv) > 1 else 'interior.leaves')
raw = path.read_bytes()
assert raw and raw.endswith(b'\n'), 'empty or unterminated certificate'
records = [line.split() for line in raw.decode('ascii').splitlines()]
assert all(len(row) == 2 for row in records), 'invalid record'
assert all(int(row[0]) in {0, 1, 2, 4, 5, 6, 7, 8} for row in records)
paths = sorted(row[1] for row in records)
assert all(p and set(p) <= {'0', '1'} for p in paths)
assert len(set(paths)) == len(paths), 'duplicate leaf'
assert all(not b.startswith(a) for a, b in zip(paths, paths[1:])), 'overlap'
depth = max(map(len, paths))
assert sum(1 << (depth - len(p)) for p in paths) == 1 << depth, 'incomplete cover'
print('PASS: prefix-free; exact Kraft sum = 1; no missing region')
print('leaves', len(paths), 'full-tree nodes', 2 * len(paths) - 1, 'maximum depth', depth)
print('methods', dict(Counter(row[0] for row in records)))
