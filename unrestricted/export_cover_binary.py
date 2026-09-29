#!/usr/bin/env python3
"""Export exact binary64 cover boxes from cover.npz for the MPFR checker."""
from __future__ import annotations
import argparse
from pathlib import Path
import struct
import numpy as np

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path, nargs='?', default=Path('cover.npz'))
    ap.add_argument('output', type=Path, nargs='?', default=Path('cover_boxes.bin'))
    args = ap.parse_args()
    data = np.load(args.source)
    val = np.asarray(data['val'], dtype='<f8', order='C')
    grad = np.asarray(data['grad'], dtype='<f8', order='C')
    classes = np.asarray(data['classes'], dtype='u1', order='C')
    assert val.ndim == grad.ndim == 2 and val.shape[1] == grad.shape[1] == 6
    assert len(classes) == len(grad)
    with args.output.open('wb') as f:
        f.write(b'SBWCOV1\0')
        f.write(struct.pack('<QQ', len(val), len(grad)))
        f.write(val.view('<u8').tobytes(order='C'))
        f.write(grad.view('<u8').tobytes(order='C'))
        f.write(classes.tobytes(order='C'))
    print(f'exported {len(val)} value boxes and {len(grad)} derivative boxes to {args.output}')

if __name__ == '__main__':
    main()
