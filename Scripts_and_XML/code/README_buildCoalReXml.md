# buildCoalReXml.R

Flexible builder for **coalRe** (BEAST 2) reassortment-network analyses. Point it
at a directory of per-segment FASTA files and it generates a ready-to-run
analysis directory (XML + LSF submission script per replicate), automatically
adapting to however many segments are present.

It is a generalization of `buildXmls.R` / `buildXmls3seg.R`: instead of
hard-coding 8 (or 3) segments, it expands/collapses the per-segment blocks of
`inference_template_wgs_cr.xml` to match the number of FASTA files found.

## Usage

```bash
Rscript buildCoalReXml.R
```

Everything is controlled by the `CONFIGURATION` block at the top of the script.
Only base R is required (no packages).

### Input requirements
- One FASTA per segment in `FASTA_DIR`.
- All segments share the **same taxon headers** (the intersection is used; a
  warning is printed if they differ, and each alignment is subset to the common
  taxa so the XML stays valid).
- Each file name ends in `_<LABEL>.fasta`; the **LABEL is the string after the
  last underscore** (e.g. `HPAI_all_North_America_HA.fasta` → `HA`). Labels
  become the segment names (`seg1`, `seg2`, … in the order given by
  `SEGMENT_ORDER`, or alphabetical by label if that is `NULL`).
- Headers are pipe-separated with the sampling date in the 3rd field
  (`name|accession|YYYY-MM-DD|…`). Adjust `DATE_SEP`/`DATE_FIELD`/`DATE_FORMAT`
  for other conventions.

### Key options
| Option | Meaning |
|---|---|
| `MODE` | `"dependent"`, `"independent"`, or `"both"` — reassortment rate coupled to Ne or not |
| `FASTA_DIR` | directory of per-segment FASTAs |
| `OUTPUT_ROOT` | output directory (holds `dependent_Ne/` and/or `independent_Ne/`) |
| `NAME` | analysis name used in output filenames |
| `N_REPS` | number of replicate folders |
| `SEGMENT_ORDER` | explicit segment ordering (labels), or `NULL` for alphabetical |
| `CLOCK_RATE` | strict-clock rate (subs/site/year) |
| `N_FINE`, `N_COARSE`, `DEEP_MAX` | Ne / reassortment-rate change-point grid |
| `SCHEDULER` | `"LSF"`, `"SLURM"`, or `"PBS"` — cluster type for the generated `sub.sh` |
| `SUB_*` | submission-script settings (walltime, memory, email, modules, beast command) |

`SUB_WALLTIME` is given LSF-style as `HH:MM`; for SLURM/PBS a `:00` is appended
automatically to make `HH:MM:SS`. `SUB_MEM` is in MB.

### `dependent` vs `independent`
- **dependent** — reassortment rate is a function of Ne (the `isNeSkyline`
  block; `SkygrowthReassortmentRatesFromSkygrowthNe` receives `logNe`). Keeps the
  log-normal prior on the reassortment rate.
- **independent** — reassortment rate is independent of Ne (the `isISkyline`
  block; no `logNe`), with a `LogDifference`/Normal smoothing prior on the
  reassortment rate.

In both cases the constant-population (`isconstant`) block is disabled.

## Output layout
```
<OUTPUT_ROOT>/
  dependent_Ne/            # (and/or independent_Ne/)
    <LABEL>.fasta          # cleaned, taxa-reconciled segment alignments
    rep1/  <NAME>.dependent.rep1.xml   sub.sh
    rep2/  ...
    rep3/  ...
```
The XMLs reference the alignments as `../<LABEL>.fasta`, so each analysis
directory is self-contained and portable.
