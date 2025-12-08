import numpy as np
import pandas as pd

from iobrpy.workflow.epic import EPIC, mRNA_cell_default


def test_epic_matches_r_fixture_for_nk_alignment():
    genes = ['g1', 'g2', 'g3']
    ref_profiles = pd.DataFrame(
        [
            [2.0, 1.0, 0.5],
            [1.0, 2.0, 0.5],
            [0.2, 0.1, 1.5],
        ],
        index=genes,
        columns=['Bcells', 'NKcells', 'Tcells'],
    )

    reference = {
        'refProfiles': ref_profiles,
        'refProfiles.var': None,
        'sigGenes': genes,
        'mRNA_cell': mRNA_cell_default,
        'var_present': False,
    }

    bulk = pd.DataFrame(
        {
            'sample1': np.array([0.775, 0.525, 0.44]),
            'sample2': np.array([1.0, 0.25, 0.8]),
        },
        index=genes,
    )

    expected = pd.read_csv('tests/fixtures/epic_r_cellfractions.csv', index_col=0)

    res = EPIC(
        bulk,
        reference,
        mRNA_cell=None,
        mRNA_cell_sub=None,
        sig_genes=genes,
        scale_exprs=False,
        with_other_cells=True,
        constrained_sum=True,
        range_based_optim=False,
    )

    cf = res['cellFractions']

    assert set(expected.columns) == set(cf.columns)
    np.testing.assert_allclose(cf.loc[expected.index, expected.columns], expected, rtol=2e-3, atol=5e-5)

    # Explicit NK regression check (ensures solver does not collapse NK to a constant column)
    assert not np.allclose(cf['NKcells'].values, cf['NKcells'].values[0])
