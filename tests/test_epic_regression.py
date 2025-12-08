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
    np.testing.assert_allclose(cf.loc['sample1', expected.columns], expected.loc['sample1'], rtol=2e-3, atol=1e-8)

    # Explicit NK regression check (ensures solver does not zero-out NK)
    assert np.isclose(cf.loc['sample1', 'NKcells'], expected.loc['sample1', 'NKcells'], rtol=2e-3)
