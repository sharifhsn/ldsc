# Python LDSC Py2 vs Py3 Notes (Minimal)

This is a minimal record of the targeted diff between `/Users/sharif/Code/ldsc/ldsc_py2` and `/Users/sharif/Code/ldsc/ldsc_py3`. Most changes are Python 3 compatibility edits and do **not** alter behavior, but a few small deltas exist.

## Potential Behavioral Deltas
1. Default LD score suffix in `parse.ldscore`.
- Py2 uses `.l2.ldscore` by default.
- Py3 prefers `.l2.ldscore.gz`.
- File: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/parse.py`.

2. Gzip handling in `read_csv`.
- Py3 adds custom gzip handling with warning interception and may return `None` or error strings instead of throwing.
- File: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/parse.py`.

3. `munge_sumstats` column alias change.
- Py2 accepts `N_CASE` as `N_CAS`.
- Py3 removes the `N_CASE` alias.
- File: `/Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py`.

4. rg denominator NaN handling.
- Py3 wraps `denom_delete_values` with `np.errstate(invalid='ignore')`, suppressing warnings for NaNs.
- File: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/regressions.py`.

5. Additional stdout noise.
- Py3 prints `cname_translation` in `munge_sumstats`.
- File: `/Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py`.

6. Py3-only helper modules.
- `ldsc_utils.py` and `ldsc_utils_local.py` exist only in py3.
- Directory: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/`.

## Non-Behavioral Edits (Python 3 Compatibility)
- `xrange` → `range`, print function syntax, relative imports.
- `np.linalg.lstsq(..., rcond=None)`.
- `bz2.BZ2File` → `bz2.open`.
- `sep='\s+'` vs `delim_whitespace=True`.
- Minor formatting and docstring corrections.

## Conclusion
No major functional differences were found. The deltas above are small and unlikely to affect typical LDSC usage, but they are documented here for completeness.
