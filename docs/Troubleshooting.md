---
title: Troubleshooting
layout: default
nav_order: 6
---

# Troubleshooting

- **Wrong input orientation**  
  Deconvolution commands expect **genes × samples**. For `deside`, `--transpose` can be helpful depending on your file.

- **Mixed separators / encoding**  
  Prefer `.csv` , `.txt` or `.tsv` consistently. Auto‑detection works in most subcommands but you can override with explicit flags where provided.

- **DeSide model missing**
  The `deside` subcommand requires pretrained model files. If you get errors like `FileNotFoundError: DeSide_model not found` , download the official model archive from:
  **https://figshare.com/articles/dataset/DeSide_model/25117862/1?file=44330255**

- **Python version for DeSide**
  The `deside` subcommand runs **ONLY on Python 3.9**. Other versions (3.8/3.10/3.11/…) are **not supported** .When invoked via the `iobrpy CLI`, it **automatically creates/uses an isolated virtual environment with pinned dependencies** so it doesn’t leak packages from your outer env. You can override the venv location with IOBRPY_DESIDE_VENV or force a clean rebuild with IOBRPY_DESIDE_REBUILD=1; the CLI wires iobrpy into that venv through a small shim and then launches the worker. 