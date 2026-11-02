# BARBIE: release-readiness review

Reviewed the supplied source snapshot on 2026-11-02. This is a source/build audit, not certification of the estimator derivation or numerical results. Line references below refer to the unchanged source. The scientific C source is preserved; proposed fixes are not silently applied.

## Blocking findings

1. **Analytic accumulator is uninitialized.** `barbie.c:87` declares `bin_Biana`; the zeroing loop around lines 201–207 initializes `bin_BS_real` and `bin_normalization` but not `bin_Biana`. Line 361 then uses `+=`. Add `bin_Biana[i][j][l][jj]=0.0` in that initialization loop. The analytic-comparison output is undefined until fixed; this finding alone does not prove the estimator numerator is wrong.
2. **Generator and estimator disagree on path and physical defaults.** `generate_field.c:98` writes under `ResFiles_N*/data_files/`, while `barbie.c:178` reads from the working directory. Generator `LL=0.56`, `ff=0.5`; estimator `LL=1`, analytic `ff=0.1`. The symlink in the README solves only the path. Choose the intended physical cell size and model parameters, apply them consistently to generator and estimator, and regenerate the field and caches. Changing a symlink does not repair units or model normalization.
3. **Cache identity omits scientific configuration.** `barbie.c:163–168` reuses `closed_rings_<NN>` based only on existence; plane-ring reuse also checks a single filename. Changes in ring definitions or field content may reuse incompatible products. Store caches by a complete configuration/input hash, or use a clean isolated run directory. Do not reuse just because one sentinel exists.
4. **CLI validation and bin bounds are insufficient.** Only argument count is checked. Width zero, reversed/indivisible intervals, or grid sizes incompatible with `Nmax=100` can yield invalid arrays/indexing or invalid estimator geometry. Require positive even NN, positive width and divisible interval; bound every bin index before accumulation. The piecewise ring-generation loop also assumes `Nr` is compatible with blocks of ten.
5. **Normalization occurs before checking the corresponding denominator.** `barbie.c:390–392` divides each multipole by its denominator, then checks only the monopole normalization. Require the specific multipole denominator to be finite and positive/nonzero as appropriate, and represent unsupported bins explicitly.

## Scientific decisions to check

- `barbie.c:355–358` sets `Num=0` whenever the measured bispectrum numerator is exactly zero. A zero measured triple product need not mean that no triangles exist. Separate geometric support from measured signal; validate any change against direct counts.
- The fixed line-of-sight, gridded field approach does not itself implement survey geometry, catalog selection, shot-noise treatment or a varying line of sight. Do not advertise these as included capabilities.
- `npoles=16` is tied to `cal_coeff`; arbitrary ell/m support is not exposed through a parameter.
- `calpow_mom` and historical variants require their own input/path/normalization review; compiling them does not select a validated reference.

## Engineering and scale

Many fopen/fread/allocation calls lack checks; a missing directory/input may crash. Binary files use native doubles and no self-description. Several large variable-length arrays are on the stack. Replace them with checked heap allocations if needed; measure peak memory/disk for the documented reference run.

The original Makefile references missing root files, invokes cleanup as part of its default target, and embeds placeholder library paths. The new root Makefile builds only available reference targets and preserves the original at `docs/legacy_makefile.txt`. Historical Slurm scripts use a different cluster's modules, paths and absent job filenames.

## Required numerical acceptance test

After fixes, generate a small but geometrically valid non-Gaussian field; compare a selected set of ring/plane triangles with explicit direct summation, including a zero-signal case and empty bins. Test all 16 output components, verify units/volume factors, and compare to the discrete model before assessing continuum differences. Add an input fixture and exact expected column definitions to the release.
