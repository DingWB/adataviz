import matplotlib

matplotlib.use(
    "Agg"
)  # non-interactive backend, must be set before any other mpl import

import numpy as np
import pandas as pd
import anndata
import pytest


@pytest.fixture
def adata():
    """Small synthetic AnnData with UMAP coords, categorical and continuous obs."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 60, 15

    X = rng.integers(0, 50, size=(n_cells, n_genes)).astype(np.float32)

    cell_types = ["TypeA", "TypeB", "TypeC", "TypeD", "TypeE"]
    cell_type_labels = [cell_types[i % len(cell_types)] for i in range(n_cells)]

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(cell_type_labels, categories=cell_types),
            "condition": pd.Categorical(
                ["ctrl" if i % 2 == 0 else "treat" for i in range(n_cells)]
            ),
            "n_genes": rng.uniform(100, 3000, size=n_cells),
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    ad = anndata.AnnData(X=X, obs=obs, var=var)
    ad.obsm["X_umap"] = rng.uniform(-10, 10, size=(n_cells, 2))

    return ad


@pytest.fixture
def adata_with_colors(adata):
    """adata with cell_type colors pre-populated in uns (avoids sc.pl.embedding call)."""
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    adata.uns["cell_type_colors"] = colors
    return adata


@pytest.fixture
def scatter_df(adata):
    """DataFrame with umap_0/umap_1 coords and a categorical hue column."""
    df = pd.DataFrame(
        {
            "umap_0": adata.obsm["X_umap"][:, 0],
            "umap_1": adata.obsm["X_umap"][:, 1],
            "cell_type": adata.obs["cell_type"].values,
            "n_genes": adata.obs["n_genes"].values,
        },
        index=adata.obs_names,
    )
    return df
