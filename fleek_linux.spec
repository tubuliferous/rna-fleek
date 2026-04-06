# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for RNA-FLEEK Linux executable"""

from pathlib import Path
from PyInstaller.utils.hooks import copy_metadata, collect_submodules

block_cipher = None

extra_datas = []
for pkg in ['scikit-learn', 'scanpy', 'anndata', 'scipy', 'numpy', 'matplotlib',
            'h5py', 'numba', 'llvmlite', 'umap-learn', 'pandas', 'packaging',
            'natsort', 'joblib', 'session-info', 'threadpoolctl', 'tqdm',
            'networkx', 'patsy', 'statsmodels', 'pynndescent', 'igraph',
            'leidenalg', 'certifi']:
    try:
        extra_datas += copy_metadata(pkg)
    except Exception:
        pass

try:
    import certifi
    extra_datas += [(certifi.where(), 'certifi')]
except ImportError:
    pass

a = Analysis(
    ['rna_fleek/server.py'],
    pathex=[],
    binaries=[],
    datas=[('rna_fleek/fleek.html', 'rna_fleek'), ('rna_fleek/cell_markers.json', 'rna_fleek')] + extra_datas,
    hiddenimports=(
        ['h5py', 'tables', 'natsort', 'packaging', 'session_info',
         'joblib', 'pandas', 'networkx', 'tqdm', 'certifi', 'ssl']
        + collect_submodules('scanpy')
        + collect_submodules('anndata')
        + collect_submodules('sklearn')
        + collect_submodules('scipy')
        + collect_submodules('matplotlib')
        + collect_submodules('umap')
        + collect_submodules('numba')
        + collect_submodules('llvmlite')
        + collect_submodules('igraph')
        + collect_submodules('leidenalg')
    ),
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'IPython', 'jupyter', 'notebook', 'pytest'],
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(pyz, a.scripts, [], exclude_binaries=True,
          name='RNA-FLEEK', debug=False, strip=False, upx=False,
          console=True)  # Show console so users see server URL and errors

coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas,
               strip=False, upx=False, name='RNA-FLEEK')
