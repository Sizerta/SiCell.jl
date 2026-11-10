#!/usr/bin/env python3
"""Generate a tiny H5AD fixture with the reference Python anndata package.

Usage: python test/fixtures/generate_python_h5ad_fixture.py
Requires: anndata, numpy, pandas, scipy
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "python_anndata_fixture.h5ad"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

# 3 observations × 4 variables. Values are intentionally non-symmetric so
# orientation mistakes are visible in the Julia regression test.
X = sparse.csr_matrix(
    np.array(
        [
            [1, 0, 3, 0],
            [0, 2, 0, 4],
            [5, 0, 0, 1],
        ],
        dtype=np.float32,
    )
)
obs = pd.DataFrame(
    {"batch": ["batch_a", "batch_a", "batch_b"], "qc_score": [0.1, 0.2, 0.3]},
    index=["py_cell_1", "py_cell_2", "py_cell_3"],
)
var = pd.DataFrame(
    {"gene_id": ["py_id_1", "py_id_2", "py_id_3", "py_id_4"]},
    index=["PY_G1", "PY_G2", "PY_G3", "PY_G4"],
)
adata = ad.AnnData(X=X, obs=obs, var=var)
adata.obsm["X_pca"] = np.array(
    [[1.0, 0.5], [2.0, 1.5], [3.0, 2.5]], dtype=np.float32
)
adata.write_h5ad(OUTPUT, compression=None)
print(f"Wrote Python anndata fixture to {OUTPUT}")
