"""Top-level entry point for `python fleek_server.py …`.

The actual implementation lives in `rna_fleek/server.py` — this is a
shim so the historical dev workflow keeps working without maintaining
a duplicate ~7000-line file. The PyInstaller specs and pip install
both already use the package layout; this file just delegates to it.
"""
from rna_fleek.server import main

if __name__ == "__main__":
    main()
