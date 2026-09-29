#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import tempfile
base = Path(__file__).resolve().parent
expected = (base / 'cover_boxes.bin').read_bytes()
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / 'cover_boxes.bin'
    subprocess.run([sys.executable, str(base / 'export_cover_binary.py'),
                    str(base / 'cover.npz'), str(out)], check=True)
    actual = out.read_bytes()
assert actual == expected, 'binary cover export differs'
print('binary cover export: byte-for-byte PASS')
