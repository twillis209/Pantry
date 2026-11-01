#!/usr/bin/env python
"""Compare two Pantry output directories (phenotyping or pheast archives).

Usage: python tests/compare_outputs.py <expected_dir> <actual_dir> --report <md>
                                       [--noise <dirA> <dirB>]

Exit status 0 iff every file matches. The optional --noise pair (two runs of the
same upstream pipeline) calibrates the tolerance for non-deterministic tensorQTL
columns: allowed |diff| = 2 x the max |diff| seen between the noise pair.
"""
import argparse
import gzip
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

RTOL = 1e-9
NOISY_COLS = ['beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'pval_perm',
              'pval_beta', 'qval', 'pval_nominal_threshold']
# Columns that depend on which variant was reported for a phenotype.
VARIANT_DEP_COLS = ['start_distance', 'end_distance', 'ma_samples', 'ma_count', 'af',
                    'slope', 'slope_se']
BYTES_SUFFIXES = ('.bed', '.txt', '.tsv', '.pos', '.list')


def _read_bytes(path: Path) -> bytes:
    if path.name.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            return f.read()
    return path.read_bytes()


def _close(a: pd.Series, b: pd.Series, atol: float = 0.0) -> pd.Series:
    """Elementwise equality (NaN == NaN); numeric within rtol/atol."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return pd.Series(np.isclose(a.astype(float), b.astype(float), rtol=RTOL,
                                    atol=atol, equal_nan=True), index=a.index)
    return (a == b) | (a.isna() & b.isna())


def _keys(df: pd.DataFrame) -> list[str]:
    return ['phenotype_id'] + (['rank'] if 'rank' in df.columns else [])


def _align(a: pd.DataFrame, b: pd.DataFrame):
    keys = _keys(a)
    a, b = a.set_index(keys), b.set_index(keys)
    if not (a.index.is_unique and b.index.is_unique):
        raise ValueError('duplicate row keys')
    common = a.index.intersection(b.index)
    n_diff = len(a.index.symmetric_difference(b.index))
    return a.loc[common], b.loc[common], n_diff


def table_noise(a: Path, b: Path, noisy_cols=NOISY_COLS) -> tuple[dict[str, float], int]:
    """Max |diff| per noisy column and number of row keys differing between two runs."""
    da, db = pd.read_csv(a, sep='\t'), pd.read_csv(b, sep='\t')
    ca, cb, n_diff = _align(da, db)
    noise = {}
    for c in noisy_cols:
        if c in ca.columns and c in cb.columns and pd.api.types.is_numeric_dtype(ca[c]):
            d = (ca[c].astype(float) - cb[c].astype(float)).abs().dropna()
            noise[c] = float(d.max()) if len(d) else 0.0
    return noise, n_diff


SUMMARY_COLS = ['pval_nominal', 'slope', 'slope_se', 'num_var']


def indep_summary(a: Path, b: Path) -> dict:
    """Summary statistics of two cis_independent tables.

    Returns {stat: (value_a, value_b, deviation, kind)}, kind in count/keys/rel.
    """
    da, db = pd.read_csv(a, sep='\t'), pd.read_csv(b, sep='\t')
    k = ['phenotype_id', 'rank']
    ka, kb = set(map(tuple, da[k].values)), set(map(tuple, db[k].values))
    out = {
        'rows': (len(da), len(db), abs(len(da) - len(db)), 'count'),
        'distinct phenotype_id': (da.phenotype_id.nunique(), db.phenotype_id.nunique(),
                                  abs(da.phenotype_id.nunique() - db.phenotype_id.nunique()), 'count'),
        'key overlap': (len(ka & kb), len(ka & kb), len(ka ^ kb), 'keys'),
    }
    m = da.merge(db, on=k, suffixes=('_a', '_b'))
    n_common = len(m)
    m = m[m['variant_id_a'] == m['variant_id_b']]
    out['variant-changed keys'] = (n_common, len(m), n_common - len(m), 'count')
    out['no same-variant keys among common keys'] = (
        n_common, len(m), int(len(m) == 0 and n_common > 0), 'flag')
    for c in SUMMARY_COLS:
        x, y = m[c + '_a'].astype(float), m[c + '_b'].astype(float)
        nan_mis = int((x.isna() != y.isna()).sum())
        out[f'{c} NaN-pattern mismatches'] = (int(x.isna().sum()), int(y.isna().sum()), nan_mis, 'count')
        ok = ~(x.isna() | y.isna())
        rel = ((x - y).abs() / x.abs().clip(lower=1e-300))[ok]
        mx = float(rel.max()) if len(rel) else 0.0
        md = float(rel.median()) if len(rel) else 0.0
        out[f'{c} max rel diff'] = (None, None, mx, 'rel')
        out[f'{c} median rel diff'] = (None, None, md, 'rel')
    return out


def indep_check(summary: dict, noise_summary: dict) -> list[tuple]:
    """Rows (stat, expected, actual, dev, noise_dev, allowed, ok)."""
    rows = []
    for name, (va, vb, dev, kind) in summary.items():
        nd = noise_summary[name][2]
        if kind == 'count':
            allowed = max(2 * nd, 1)
        elif kind == 'keys':
            allowed = 2 * nd
        elif kind == 'rel':
            allowed = max(2 * nd, 1e-6)
        elif kind == 'flag':
            allowed = 0
        else:
            rows.append((name, va, vb, dev, nd, None, True))
            continue
        rows.append((name, va, vb, dev, nd, allowed, dev <= allowed))
    return rows


def _compare_tables(a: Path, b: Path, exact_cols, noisy_cols, noise, rowset_noise=0):
    mism, notes = [], []
    da, db = pd.read_csv(a, sep='\t'), pd.read_csv(b, sep='\t')
    if list(da.columns) != list(db.columns):
        return [f'columns differ: {list(da.columns)} vs {list(db.columns)}'], notes
    try:
        ca, cb, n_diff = _align(da, db)
    except ValueError as e:
        return [str(e)], notes
    if n_diff:
        if 0 < n_diff <= 2 * rowset_noise:
            notes.append(f'noise-level row-set difference: {n_diff} differing keys '
                         f'(noise pair: {rowset_noise})')
        else:
            mism.append(f'row set differs: {n_diff} differing keys '
                        f'(allowed up to {2 * rowset_noise})')
    tie = pd.Series(False, index=ca.index)
    if 'variant_id' in ca.columns:
        vdiff = ~_close(ca['variant_id'], cb['variant_id'])
        if 'pval_nominal' in ca.columns:
            same_p = _close(ca['pval_nominal'], cb['pval_nominal'])
        else:
            same_p = pd.Series(False, index=ca.index)
        tie = vdiff & same_p
        if (vdiff & ~same_p).any():
            ex = ca.index[(vdiff & ~same_p).values][:3].tolist()
            mism.append(f'variant_id differs with different pval_nominal in '
                        f'{int((vdiff & ~same_p).sum())} rows, e.g. {ex}')
        if tie.any():
            notes.append(f'tie-break? variant_id differs with equal pval_nominal in '
                         f'{int(tie.sum())} rows')
    for c in da.columns:
        if c in _keys(da):
            continue
        if c == 'variant_id':
            continue
        is_noisy = c in noisy_cols and pd.api.types.is_numeric_dtype(ca[c])
        if is_noisy:
            tol = 2 * noise.get(c, 0.0)
            bad = ~_close(ca[c], cb[c], atol=tol)
            what = f'exceeds allowed noise {tol:.3g}'
        else:
            bad = ~_close(ca[c], cb[c])
            if c in exact_cols and c in VARIANT_DEP_COLS:
                bad &= ~tie
            what = 'differs'
        if bad.any():
            ex = ca.index[bad.values][0]
            mism.append(f'{c} {what} in {int(bad.sum())} rows, e.g. {ex}: '
                        f'{ca.loc[ex, c]} vs {cb.loc[ex, c]}')
    return mism, notes


def compare_tables(a: Path, b: Path, exact_cols, noisy_cols, noise, rowset_noise=0) -> list[str]:
    """Return a list of mismatch descriptions (empty if the tables agree)."""
    return _compare_tables(a, b, exact_cols, noisy_cols, noise, rowset_noise)[0]


def _tensorqtl_cols(path: Path):
    cols = pd.read_csv(path, sep='\t', nrows=0).columns
    return [c for c in cols if c not in NOISY_COLS], [c for c in NOISY_COLS if c in cols]


def _compare_rdat(pairs: list[tuple[Path, Path]]) -> list[str]:
    if not pairs:
        return []
    mism = []
    script = ('a<-new.env();b<-new.env();load(A,a);load(B,b);'
              'stopifnot(isTRUE(all.equal(as.list(a),as.list(b),tolerance=1e-8)))')
    for x, y in pairs:
        try:
            r = subprocess.run(['Rscript', '-e', f'A<-"{x}";B<-"{y}";{script}'],
                               capture_output=True, text=True)
        except OSError:
            return ['RDat not compared: Rscript unavailable']
        if r.returncode != 0:
            mism.append(f'{x.name}: RDat differs ({r.stderr.strip().splitlines()[-1:]})')
    return mism


def _compare_tar(a: Path, b: Path):
    mism = []
    with tempfile.TemporaryDirectory() as td:
        roots = []
        for side, p in (('a', a), ('b', b)):
            root = Path(td) / side
            root.mkdir()
            with tarfile.open(p) as t:
                t.extractall(root, filter='data')
            roots.append(root)
        fa = {p.relative_to(roots[0]) for p in roots[0].rglob('*') if p.is_file()}
        fb = {p.relative_to(roots[1]) for p in roots[1].rglob('*') if p.is_file()}
        for r in sorted(fa ^ fb):
            mism.append(f'{r}: present on one side only')
        rdat = []
        for r in sorted(fa & fb):
            n = r.name
            if n.endswith('.profile.err'):
                continue
            if n.endswith('.wgt.RDat'):
                rdat.append((roots[0] / r, roots[1] / r))
            elif (roots[0] / r).read_bytes() != (roots[1] / r).read_bytes():
                mism.append(f'{r}: contents differ')
        mism += _compare_rdat(rdat)
    return mism


def _bed_numeric(a: Path, b: Path):
    """Numeric BED comparison. Returns (errors, max |diff| over non-NaN value cells)."""
    la = _read_bytes(a).decode().splitlines()
    lb = _read_bytes(b).decode().splitlines()
    if la[:1] != lb[:1]:
        return ['header differs'], None
    if len(la) != len(lb):
        return [f'row count differs ({len(la) - 1} vs {len(lb) - 1})'], None
    maxd, nan_bad, inf_bad, id_bad, bad_val = 0.0, 0, 0, [], []
    for i, (x, y) in enumerate(zip(la[1:], lb[1:]), start=1):
        fx, fy = x.split('\t'), y.split('\t')
        if len(fx) != len(fy) or fx[:4] != fy[:4]:
            id_bad.append(i)
            continue
        try:
            vx = np.array([float(v) if v not in ('NA', 'nan', 'NaN', '') else np.nan for v in fx[4:]])
            vy = np.array([float(v) if v not in ('NA', 'nan', 'NaN', '') else np.nan for v in fy[4:]])
        except ValueError:
            bad_val.append(i)
            continue
        if (np.isnan(vx) != np.isnan(vy)).any():
            nan_bad += 1
            continue
        ok = ~np.isnan(vx)
        vx, vy = vx[ok], vy[ok]
        same = vx == vy  # equal values, including matching infinities
        if (np.isinf(vx) | np.isinf(vy))[~same].any():
            inf_bad += 1
            continue
        if (~same).any():
            maxd = max(maxd, float(np.abs(vx[~same] - vy[~same]).max()))
    errs = []
    if id_bad:
        errs.append(f'first 4 columns (or column count) differ in {len(id_bad)} rows, e.g. line {id_bad[0]}')
    if nan_bad:
        errs.append(f'NaN pattern differs in {nan_bad} rows')
    if inf_bad:
        errs.append(f'infinite values differ in {inf_bad} rows')
    if bad_val:
        errs.append(f'unparseable value in row {bad_val[0]} ({len(bad_val)} rows)')
    if not errs and np.isnan(maxd):
        errs.append('max |diff| is NaN')
    return errs, (None if errs else maxd)


def _compare_hsq(a: Path, b: Path):
    la, lb = a.read_text().splitlines(), b.read_text().splitlines()
    if la[:1] != lb[:1]:
        return ['header differs']
    if sorted(la[1:]) != sorted(lb[1:]):
        return [f'rows differ ({len(set(la[1:]) ^ set(lb[1:]))} differing lines)']
    return []


def _files(d: Path) -> set[str]:
    return {str(p.relative_to(d)) for p in d.rglob('*') if p.is_file()}


def compare_dirs(expected: Path, actual: Path, noise_dirs=None, summaries=None):
    """Compare two dirs. Returns (mismatches: {relpath: [msgs]} flattened to list, notes)."""
    expected, actual = Path(expected), Path(actual)
    bad = [f'{d} is not an existing directory' for d in
           (expected, actual, *(noise_dirs or ())) if not Path(d).is_dir()]
    if bad:
        return bad, []
    fe, fa = _files(expected), _files(actual)
    mism, notes = [], []
    if not fe:
        mism.append(f'expected directory {expected} contains no files')
    for r in sorted(fe - fa):
        mism.append(f'{r}: missing from actual')
    for r in sorted(fa - fe):
        mism.append(f'{r}: missing from expected (extra in actual)')
    for r in sorted(fe & fa):
        e, a = expected / r, actual / r
        name = r
        if name.endswith('.tbi'):
            continue  # presence checked above; content covered by the BED it indexes
        msgs, nts = [], []
        if name.endswith(('.cis_qtl.txt.gz', '.cis_independent_qtl.txt.gz')):
            exact, noisy = _tensorqtl_cols(e)
            noise, rs = {}, 0
            if noise_dirs is not None:
                na, nb = Path(noise_dirs[0]) / r, Path(noise_dirs[1]) / r
                if na.exists() and nb.exists():
                    noise, rs = table_noise(na, nb)
                else:
                    mism.append(f'{r}: noise counterpart missing ({na if not na.exists() else nb})')
                    continue
            if 'independent' not in name:
                rs = 0
            if 'independent' in name:
                summ = indep_summary(e, a)
                nsumm = None
                if noise_dirs is not None:
                    nsumm = indep_summary(na, nb)
                    chk = indep_check(summ, nsumm)
                    msgs = [f'summary {c[0]}: deviation {c[3]:.3g} > allowed {c[5]:.3g} '
                            f'(noise {c[4]:.3g})' for c in chk if not c[6]]
                else:
                    chk = [(n, v[0], v[1], v[2], None, None, True) for n, v in summ.items()]
                    msgs, nts = _compare_tables(e, a, exact, noisy, noise, rs)
                if summaries is not None:
                    summaries.append((r, chk))
            else:
                msgs, nts = _compare_tables(e, a, exact, noisy, noise, rs)
        elif name.endswith('.tar.bz2'):
            msgs = _compare_tar(e, a)
        elif name.endswith('.hsq.tsv'):
            msgs = _compare_hsq(e, a)
        elif _read_bytes(e) != _read_bytes(a):
            msgs = ['contents differ']
            if noise_dirs is not None and name.endswith(('.bed', '.bed.gz')):
                na, nb = Path(noise_dirs[0]) / r, Path(noise_dirs[1]) / r
                if not (na.exists() and nb.exists()):
                    msgs = [f'noise counterpart missing ({na if not na.exists() else nb})']
                else:
                    nerr, nmax = _bed_numeric(na, nb)
                    if nerr:
                        msgs = [f'noise pair not numerically comparable: {"; ".join(nerr)}']
                    else:
                        errs, dmax = _bed_numeric(e, a)
                        allowed = 2 * nmax
                        if errs:
                            msgs = errs
                        elif dmax > allowed:
                            msgs = [f'max |diff| {dmax:.4g} > allowed {allowed:.4g} (noise pair max {nmax:.4g})']
                        else:
                            msgs = []
                            nts = [f'numeric BED: max |diff| {dmax:.4g}, allowed {allowed:.4g}']
        mism += [f'{r}: {m}' for m in msgs]
        notes += [f'{r}: {m}' for m in nts]
    return mism, notes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('expected', type=Path)
    ap.add_argument('actual', type=Path)
    ap.add_argument('--report', type=Path, required=True)
    ap.add_argument('--noise', nargs=2, type=Path, metavar=('DIR_A', 'DIR_B'))
    args = ap.parse_args(argv)
    summaries = []
    mism, notes = compare_dirs(args.expected, args.actual, args.noise, summaries)
    lines = [f'# Comparison: {args.expected} vs {args.actual}', '',
             f'Noise pair: {" , ".join(map(str, args.noise)) if args.noise else "none"}', '',
             ('With --noise, cis_independent files are compared ONLY on the summary statistics below '
              '(no row-level comparison).' if args.noise else
              'Without --noise, cis_independent files are compared row by row (exact).'), '',
             f'## Mismatches ({len(mism)})', ''] + [f'- {m}' for m in mism] + \
            ['', f'## Notes ({len(notes)})', ''] + [f'- {n}' for n in notes] + \
            ['', '## cis_independent summaries (expected | actual | deviation | noise deviation | allowed | ok)', ''] + \
            sum([[f'### {r}', '', '| statistic | expected | actual | deviation (rel rows: rel diff vs expected) | noise dev (rel rows: noise-pair rel diff) | allowed | ok |',
                  '|---|---|---|---|---|---|---|'] +
                 [f'| {c[0]} | {"-" if c[1] is None else format(c[1], ".4g")} | {"-" if c[2] is None else format(c[2], ".4g")} | {c[3]:.3g} | '
                  f'{"-" if c[4] is None else format(c[4], ".3g")} | '
                  f'{"-" if c[5] is None else format(c[5], ".3g")} | {"yes" if c[6] else "NO"} |' for c in chk] + ['']
                 for r, chk in summaries], []) + \
            ['', f'**Result: {"MATCH" if not mism else "MISMATCH"}**', '']
    args.report.write_text('\n'.join(lines))
    print('\n'.join(lines))
    return 0 if not mism else 1


if __name__ == '__main__':
    sys.exit(main())
