import gzip
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import compare_outputs as co  # noqa: E402

COLS = ['phenotype_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df',
        'variant_id', 'start_distance', 'end_distance', 'ma_samples', 'ma_count', 'af',
        'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta', 'qval',
        'pval_nominal_threshold']


def cis_df(**over):
    df = pd.DataFrame({
        'phenotype_id': ['g1', 'g2'], 'num_var': [10, 20],
        'beta_shape1': [1.0, 1.1], 'beta_shape2': [300.0, 310.0],
        'true_df': [350.0, 360.0], 'pval_true_df': [0.01, 0.02],
        'variant_id': ['chr1_1_A_C_b38', 'chr1_2_A_C_b38'],
        'start_distance': [5, 6], 'end_distance': [5, 6], 'ma_samples': [20, 21],
        'ma_count': [22, 23], 'af': [0.1, 0.2], 'pval_nominal': [0.001, 0.002],
        'slope': [0.5, -0.4], 'slope_se': [0.1, 0.2], 'pval_perm': [0.5, 0.6],
        'pval_beta': [0.5, 0.6], 'qval': [0.3, 0.4],
        'pval_nominal_threshold': [0.0001, 0.0002]})
    for k, v in over.items():
        df[k] = v
    return df[COLS]


def write_gz(path, df):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt') as f:
        df.to_csv(f, sep='\t', index=False)


def write_bed(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt') as f:
        f.write(text)


def make_pair(tmp_path, name, a, b):
    for side, df in (('exp', a), ('act', b)):
        write_gz(tmp_path / side / 'qtl' / name, df)
    return tmp_path / 'exp', tmp_path / 'act'


def run(tmp_path, a, b, noise=None, name='x.cis_qtl.txt.gz'):
    exp, act = make_pair(tmp_path, name, a, b)
    nd = None
    if noise is not None:
        na, nb = noise
        write_gz(tmp_path / 'na' / 'qtl' / name, na)
        write_gz(tmp_path / 'nb' / 'qtl' / name, nb)
        nd = (tmp_path / 'na', tmp_path / 'nb')
    return co.compare_dirs(exp, act, nd)


def test_identical_bed_passes(tmp_path):
    for s in ('exp', 'act'):
        write_bed(tmp_path / s / 'e.bed.gz', '#chr\tstart\n1\t2\n')
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert mism == []


def test_bed_value_change_fails(tmp_path):
    write_bed(tmp_path / 'exp' / 'e.bed.gz', '#chr\tstart\n1\t2\n')
    write_bed(tmp_path / 'act' / 'e.bed.gz', '#chr\tstart\n1\t3\n')
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert len(mism) == 1 and 'e.bed.gz' in mism[0]


def test_noisy_column_within_noise_passes(tmp_path):
    base = cis_df()
    noise = (base, cis_df(pval_perm=[0.51, 0.6]))
    mism, _ = run(tmp_path, base, cis_df(pval_perm=[0.515, 0.6]), noise)
    assert mism == []


def test_noisy_column_outside_noise_fails(tmp_path):
    base = cis_df()
    noise = (base, cis_df(pval_perm=[0.51, 0.6]))
    mism, _ = run(tmp_path, base, cis_df(pval_perm=[0.53, 0.6]), noise)
    assert any('pval_perm' in m for m in mism)


def test_noisy_column_without_noise_pair_fails(tmp_path):
    mism, _ = run(tmp_path, cis_df(), cis_df(qval=[0.31, 0.4]))
    assert any('qval' in m for m in mism)


def test_exact_column_change_fails(tmp_path):
    mism, _ = run(tmp_path, cis_df(), cis_df(slope=[0.5, -0.41]))
    assert any('slope' in m for m in mism)


def test_tie_break_reported_not_failed(tmp_path):
    other = cis_df(variant_id=['chr1_9_A_C_b38', 'chr1_2_A_C_b38'],
                   start_distance=[8, 6], slope=[0.9, -0.4])
    mism, notes = run(tmp_path, cis_df(), other)
    assert mism == []
    assert any('tie-break?' in n for n in notes)


def test_variant_change_with_different_pval_fails(tmp_path):
    other = cis_df(variant_id=['chr1_9_A_C_b38', 'chr1_2_A_C_b38'],
                   pval_nominal=[0.003, 0.002])
    mism, _ = run(tmp_path, cis_df(), other)
    assert any('variant_id' in m for m in mism)


def test_missing_file_fails(tmp_path):
    write_bed(tmp_path / 'exp' / 'e.bed.gz', 'a\n')
    (tmp_path / 'act').mkdir()
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert len(mism) == 1 and 'missing' in mism[0]


def indep(rows):
    df = cis_df().iloc[[0, 1]].copy()
    df = pd.concat([df.iloc[[0]]] * len(rows)).reset_index(drop=True)
    df['phenotype_id'] = [r[0] for r in rows]
    df['rank'] = [r[1] for r in rows]
    return df.drop(columns=['qval', 'pval_nominal_threshold'])


def test_cis_independent_rowset_diff_noise_level(tmp_path):
    name = 'x.cis_independent_qtl.txt.gz'
    a = indep([('g1', 1), ('g2', 1)])
    b = indep([('g1', 1), ('g2', 1), ('g2', 2)])
    noise = (a, b)  # one differing key
    mism, notes = run(tmp_path, a, b, noise, name)
    assert mism == []
    assert any('noise-level row-set difference' in n for n in notes)


def test_cis_independent_rowset_diff_too_big_fails(tmp_path):
    name = 'x.cis_independent_qtl.txt.gz'
    a = indep([('g1', 1)])
    b = indep([('g1', 1), ('g2', 1), ('g2', 2), ('g3', 1)])
    noise = (a, indep([('g1', 1), ('g2', 1)]))  # 1 differing key; actual has 3
    mism, _ = run(tmp_path, a, b, noise, name)
    assert any('row set' in m for m in mism)


def test_cis_independent_rowset_diff_without_noise_fails(tmp_path):
    name = 'x.cis_independent_qtl.txt.gz'
    a = indep([('g1', 1)])
    b = indep([('g1', 1), ('g1', 2)])
    mism, _ = run(tmp_path, a, b, None, name)
    assert any('row set' in m for m in mism)


def test_cis_independent_rowset_diff_noise_pair_identical_fails(tmp_path):
    name = 'x.cis_independent_qtl.txt.gz'
    a = indep([('g1', 1)])
    b = indep([('g1', 1), ('g1', 2)])
    mism, _ = run(tmp_path, a, a, None, name)  # sanity: identical passes
    assert mism == []
    mism, _ = run(tmp_path, a, b, (a, a), name)
    assert any('row set' in m for m in mism)


def test_heritability_row_order_ignored(tmp_path):
    for s, txt in (('exp', 'phenotype_id\thsq\na\t1\nb\t2\n'), ('act', 'phenotype_id\thsq\nb\t2\na\t1\n')):
        p = tmp_path / s / 'heritability' / 'x.hsq.tsv'
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(txt)
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert mism == []


def test_tbi_follows_presence(tmp_path):
    write_bed(tmp_path / 'exp' / 'e.bed.gz', 'a\n')
    write_bed(tmp_path / 'act' / 'e.bed.gz', 'a\n')
    (tmp_path / 'exp' / 'e.bed.gz.tbi').write_bytes(b'x')
    (tmp_path / 'act' / 'e.bed.gz.tbi').write_bytes(b'different compressed bytes')
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert mism == []
    (tmp_path / 'act' / 'e.bed.gz.tbi').unlink()
    mism, _ = co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')
    assert len(mism) == 1


def test_twas_tar_profile_err_ignored_and_list_compared(tmp_path):
    import tarfile
    import io

    def mk(side, lst, err):
        p = tmp_path / side / 'twas' / 'x.tar.bz2'
        p.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(p, 'w:bz2') as t:
            for n, d in (('x.list', lst), ('x.profile.err', err)):
                ti = tarfile.TarInfo(n)
                ti.size = len(d)
                t.addfile(ti, io.BytesIO(d))
    mk('exp', b'a\n', b'1s')
    mk('act', b'a\n', b'2s')
    assert co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')[0] == []
    mk2 = tmp_path / 'act' / 'twas' / 'x.tar.bz2'
    mk2.unlink()
    mk('act', b'b\n', b'2s')
    assert len(co.compare_dirs(tmp_path / 'exp', tmp_path / 'act')[0]) == 1
