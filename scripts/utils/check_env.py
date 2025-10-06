#!/usr/bin/env python3
import os, sys, yaml

cfg_path = 'config/paths.yaml'
if not os.path.exists(cfg_path):
    print('Missing config/paths.yaml (copy paths.example.yaml).', file=sys.stderr); sys.exit(1)
cfg = yaml.safe_load(open(cfg_path))
ok = True
for k in ['slicer','atlases','io']:
    if k not in cfg: print(f'Missing section: {k}', file=sys.stderr); ok=False
print('OK' if ok else 'NOT OK')
